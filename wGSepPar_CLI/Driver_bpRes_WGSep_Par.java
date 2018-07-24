package wGSepPar_CLI;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class Driver_bpRes_WGSep_Par {
	/** Window size for calculating cMBF (in bp) */
	private static Integer windowbpSize;
	/** Read resolution / width of bin */
	private static int intervalSize;
	/** Median multiple */
	private static double medianMult;
	/** Minimum read count (in place of zero) */
	private static double minRC;
//	private static boolean toIntegrate;
	
	private static String inBed;
	
	private static String baseName;
	private static String dirName;
	private static String outBaseName;
	private static String spacer;
	
	/**
	 * Main method
	 *@param args - interval size, window size, medianMult, minRC, {input .bed, output .bedGraph} filenames
	 * !! Window size must be a multiple of interval size !!
	 * !! Ignores any lines that do not parse as a data line (headers, etc) !!
	 */
	public static void main(String[] args) throws Exception{
		
		parseOptions(args);
		
		//split chromosome files
		ArrayList<File> chromFiles = null;
		try {
			chromFiles = splitChromFiles(inBed);
		} catch (FileNotFoundException e) {
			System.out.println(e.getMessage());	
			System.exit(1);
		}
		
		//parallel process chromosomes 
		int nThreads = Runtime.getRuntime().availableProcessors();
		ExecutorService executor = Executors.newFixedThreadPool(nThreads);
		
		ArrayList<Future<File>> futureList = new ArrayList<Future<File>>();
		ArrayList<File> outFileList = new ArrayList<File>();
		
		for (File f : chromFiles) {
			Callable<File> worker = new ProcessChromFile(dirName, f);
			Future<File> submit = executor.submit(worker);
			futureList.add(submit);
		}
		
		for (Future<File> future : futureList) {
            try {
            	File f = future.get();
                System.out.println(f.getPath() + " completed");
            	outFileList.add(f);
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
		
		executor.shutdown();
		executor.awaitTermination(12, TimeUnit.HOURS);
		
		try {
			mergeWholeChrFile(outFileList);
		} catch (FileNotFoundException e) {
			System.err.println(e.getMessage());
		}
		
		System.out.println("completed");
	}
	
	public static void parseOptions(String[] args) {
		Options options = new Options();
		
		Option infile = new Option("i", "input", true, "[req] input file path, must be a bed file");
		infile.setRequired(true);
		options.addOption(infile);
		
		Option outfileName = new Option("o", "output", true, "[opt] output file base name");
		outfileName.setRequired(false);
		options.addOption(outfileName);
		
		Option winSize = new Option("w", "windowSize", true, "[req] window size for calculating cMBF, must be a multiple of interval size");
		winSize.setRequired(true);
		options.addOption(winSize);
		
		Option medMult = new Option("m", "medMult", true, "[opt] median multiple, cannot be 0");
		medMult.setRequired(false);
		options.addOption(medMult);
		
		Option defZero = new Option("z", "defaultZero", true, "[opt] default number to replace zero");
		defZero.setRequired(false);
		options.addOption(defZero);
		
		Option help = new Option("h", "help", false, "");
		options.addOption(help);
			
		CommandLineParser parser = new DefaultParser();
		HelpFormatter formatter = new HelpFormatter();
		CommandLine cmd = null;
		
		String header = "Normalize NGS genome read count data using cMBFs \n\n"; /////////////////////////////
		String footer = "\nPlease report issues at http://example.com/issues"; ///////////////////////////////
		
		try {
			cmd = parser.parse(options, args);
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			formatter.printHelp("Driver", header, options, footer, false);
			
			System.exit(1);
		}

		if (cmd.hasOption("h")) {
			formatter.printHelp("Driver", header, options, footer, false);
		}
		
		//default options
		outBaseName = "out";
		windowbpSize = 10000;
		medianMult = 1;
		minRC = 0.5;
		
		//file options
		inBed = cmd.getOptionValue("i");
		if (cmd.hasOption("o")) {
			outBaseName = cmd.getOptionValue("o");
		}
		if (!inBed.endsWith(".bed")) {
			throw new IllegalArgumentException("Incorrect filetypes");
		}
		
		//other parameters
		if (cmd.hasOption("w")) {
			windowbpSize = Integer.parseInt(cmd.getOptionValue("w"));
		}
		if (cmd.hasOption("m")) {
			double mm = Double.parseDouble(cmd.getOptionValue("m"));
			if (mm == 0) {
				System.err.println("Median multiple cannot be zero.");
				System.exit(0);
			}
			medianMult = mm;
		}
		if (cmd.hasOption("z")) {
			double mrc = Double.parseDouble(cmd.getOptionValue("z"));
			if (mrc == 0) {
				System.err.println("Default zero cannot be zero.");
				System.exit(0);
			}
			minRC = mrc;
		}
		
		//set static parameters
		IntStats_WGSep.setMinRC(minRC);
		ProcessChromFile.setWinSize(windowbpSize);
		ProcessChromFile.setIntSize(intervalSize);
		ProcessChromFile.setMedMult(medianMult);
		ProcessChromFile.setOutBaseName(outBaseName);
	}
	
	public static ArrayList<File> splitChromFiles(String infile) throws FileNotFoundException {
		ArrayList<File> chromFiles = new ArrayList<File>();
		Scanner sc;
		try {
			sc = new Scanner(new File(infile));
		} catch (FileNotFoundException e) {
			throw new FileNotFoundException("File not found: " + infile);
		}
		
		spacer = "/";
		baseName = infile.substring(0, infile.length() - 4);
		dirName = baseName + "_ChromFiles";
		
		File dir = new File(dirName);
		boolean folderMade = dir.mkdir();
		
		if (!dir.exists() && !folderMade) {
			spacer = "_";
		}
		
		if (!sc.hasNextLine()) {
			sc.close();
			return null;
		}
		
		String line = null;
		String chrom = null;
		while (chrom == null && sc.hasNextLine()) {
			line = sc.nextLine();
			chrom = parseLineChrom(line);
		}
		
		File cFile = new File(dirName + spacer + chrom + ".bed");
		chromFiles.add(cFile);
		PrintWriter pw = new PrintWriter(cFile);
		pw.println(line);
		
		while(sc.hasNextLine()) { //separate files
			String nChrom = null;
			while (nChrom == null && sc.hasNextLine()) {
				line = sc.nextLine();
				nChrom = parseLineChrom(line);
			}
			if (!chrom.equals(nChrom)) {
				pw.close();
				chrom = nChrom;
				cFile = new File(dirName + spacer + nChrom + ".bed");
				chromFiles.add(cFile);
				pw = new PrintWriter(cFile);
			}
			pw.println(line);
		}
		
		pw.close();
		sc.close();
		
		return chromFiles;
	}
	
	/** checks if the line is able to be parsed - if not, the program will ignore the line entirely */
	public static String parseLineChrom(String line) {
		String[] ra = line.split("\\s+");

		if (ra[0].startsWith("chr")) {
			return ra[0];
		}
		
		try {
			Integer.parseInt(ra[1]);
			Integer.parseInt(ra[2]);
			Integer.parseInt(ra[3]);
		} catch (NumberFormatException e) {
			return null;
		} catch (IndexOutOfBoundsException e) {
			return null;
		}
		
		return null;
	}
	
	public static void mergeWholeChrFile(ArrayList<File> outChrFiles) throws FileNotFoundException {
		String wChrFileName = dirName.substring(0, dirName.lastIndexOf("_ChromFiles")) + "_out" + spacer + outBaseName + "_allChr.bedGraph";
		File wChrFile = new File (wChrFileName);
		PrintWriter wholePW = null;
		try {
			wholePW = new PrintWriter(wChrFile);
		} catch (FileNotFoundException e) {
			throw new FileNotFoundException("File not found: " + wChrFileName);
		}
		
		wholePW.println("track type=bedGraph name=\"" + wChrFileName + "\"" + " description=\"" + wChrFileName + "\" "
				+ "visibility=full autoScale=Off alwaysZero=On maxHeightPixels=128:30:11 viewLimits=0:1"); //header
		
		for(File cFile : outChrFiles) {
			Scanner sc = new Scanner(cFile);
			sc.nextLine(); //skip header
			while(sc.hasNextLine()) {
				wholePW.println(sc.nextLine());
			}
			sc.close();
		}
		
		wholePW.close();
		System.out.println("Merging completed");
	}
	
}
