

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
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

public class GMIM_Pipeline {

	/** Window size for calculating cMBF (in bp) */
	private static Integer windowbpSize;
	/** Read resolution / width of bin */
	private static int intervalSize;
	/** Median multiple */
	private static double medianMult;
	/** Minimum read count (in place of zero) */
	private static double minRC;

	private static String[] inBedFileNames;
	private static String[] outBaseNames;
	private static String integBaseName;

	public static void main(String[] args) throws Exception {
		parseMultOptions(args);

		//available processors / samples for each
		int nThreads = Runtime.getRuntime().availableProcessors() / inBedFileNames.length;
		ExecutorService executor = Executors.newFixedThreadPool(nThreads);
		
		ArrayList<Future<ArrayList<File>>> futureList = new ArrayList<Future<ArrayList<File>>>();
		/** (# samples) of ArrayLists of the output files*/
		ArrayList<ArrayList<File>> allOutFiles = new ArrayList<ArrayList<File>>();

		for (int i = 0; i < inBedFileNames.length; i++) {
			String inBed = inBedFileNames[i];
			String outBase = outBaseNames[i];
			
			Callable<ArrayList<File>> worker = new Driver(inBed, outBase);
			Future<ArrayList<File>> submit = executor.submit(worker);
			futureList.add(submit);
		}

		int i = 0;
		for (Future<ArrayList<File>> future : futureList) {
	        try {
	        	ArrayList<File> f = future.get();
	            System.out.println(inBedFileNames[i] + " completed\n");
	        	allOutFiles.add(f);
	        	i++;
	        } catch (InterruptedException e) {
	            e.printStackTrace();
	        } catch (ExecutionException e) {
	            e.printStackTrace();
	        }
	    }
		
		executor.shutdown();
		executor.awaitTermination(12, TimeUnit.HOURS);
		
		System.out.println("\ncMBF calculations completed\n***\n");
		
		//integrate
		if (allOutFiles.size() > 1) {
			integrateAll(allOutFiles);
			System.out.println("\nIntegration completed\n***\n");
		}
	}

	public static void parseMultOptions(String[] args) {
		Options options = new Options();

		Option infile = new Option("i", "input", true, "[req] input file paths, must be bed files");
		infile.setRequired(true);
		infile.setArgs(Option.UNLIMITED_VALUES);
		options.addOption(infile);

		Option outfilename = new Option("o", "output", true, "[opt] output file base names, must have one for each input file");
		outfilename.setRequired(false);
		outfilename.setArgs(Option.UNLIMITED_VALUES);
		options.addOption(outfilename);

		Option intfilename = new Option("int", "integrate", true, "[req] integration folder/file base name");
		intfilename.setRequired(true);
		intfilename.setArgs(1);
		options.addOption(intfilename);

		Option winSize = new Option("w", "windowSize", true, "[req] window size for calculating cMBF, must be a multiple of interval size");
		winSize.setRequired(true);
		winSize.setArgs(1);
		options.addOption(winSize);

		Option medMult = new Option("m", "medMult", true, "[opt] median multiple, cannot be 0");
		medMult.setRequired(false);
		medMult.setArgs(1);
		options.addOption(medMult);

		Option defZero = new Option("z", "defaultZero", true, "[opt] default number to replace zero");
		defZero.setRequired(false);
		defZero.setArgs(1);
		options.addOption(defZero);

		Option help = new Option("h", "help", false, "");
		options.addOption(help);

		CommandLineParser parser = new DefaultParser();
		HelpFormatter formatter = new HelpFormatter();
		CommandLine cmd = null;

		String header = "Normalize and integrate NGS genome read count data using cMBFs \n\n"; 
		String footer = "\nPlease report issues at http://example.com/issues"; 

		try {
			cmd = parser.parse(options, args);
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			formatter.printHelp("GMIM", header, options, footer, false);
			System.exit(1);
		}

		if (cmd.hasOption("h")) {
			formatter.printHelp("GMIM", header, options, footer, false);
		}

		//default options
		windowbpSize = 10000;
		medianMult = 1;
		minRC = 0.5;

		//file options
		inBedFileNames = cmd.getOptionValues("i");
		for (String inBed : inBedFileNames) {
			if (!inBed.endsWith(".bed")) {
				throw new IllegalArgumentException("Must be a .bed");
			}
		}
		if (cmd.hasOption("o")) {
			outBaseNames = cmd.getOptionValues("o");
			if (inBedFileNames.length != outBaseNames.length) {
				throw new IllegalArgumentException("Number of input files and output base names do not match");
			}
			for (String outName : outBaseNames) {
				if (outName.contains("/")) throw new IllegalArgumentException("Cannot use a path as an output base name");
			}
		}
		else {
			outBaseNames = new String[inBedFileNames.length];
			Arrays.fill(outBaseNames, "out");
		}
		if (cmd.hasOption("int")) {
			integBaseName = cmd.getOptionValue("int");
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
		Driver.windowbpSize = windowbpSize;
		Driver.intervalSize = intervalSize;
		Driver.medianMult = medianMult;
		Driver.minRC = minRC;

		IntStats_WGSep.setMinRC(minRC);
		ProcessChromFile.setWinSize(windowbpSize);
		ProcessChromFile.setIntSize(intervalSize);
		ProcessChromFile.setMedMult(medianMult);
	}

	public static void integrateAll(ArrayList<ArrayList<File>> allOutFiles) throws InterruptedException {
		//make directory for integration files
		File integDir = new File(integBaseName);
		boolean folderMade = integDir.mkdirs();

		String spacer = "/";

		if (!integDir.exists() && !folderMade) {
			spacer = "_";
		}

		//for each file in the first ArrayList, check that there are corresponding files
		//in the other ArrayLists and integrate them in parallel

		ArrayList<File> ra0 = allOutFiles.get(0);

		int nThreads = Runtime.getRuntime().availableProcessors();
		ExecutorService executor = Executors.newFixedThreadPool(nThreads);

		ArrayList<Future<File>> futureList = new ArrayList<Future<File>>();
		ArrayList<File> outFileList = new ArrayList<File>();

		for (File f : ra0) {
			//check each file in first ArrayList
			String chrName = parseBGChrom(f);

			ArrayList<File> matchChrFiles = new ArrayList<File>();
			matchChrFiles.add(f);

			//search other ArrayLists for same chromosome file
			for (int i = 1; i < allOutFiles.size(); i++) { //for each other ArrayList
				ArrayList<File> ra = allOutFiles.get(i);

				for (File fother : ra) {
					String chrOther = parseBGChrom(fother);
					if (chrName.equals(chrOther)) {
						matchChrFiles.add(fother);
						break;
					}
				}
			} 

			//found all available matching chromosome arrays
			if (matchChrFiles.size() == allOutFiles.size()) { //all samples have a chromosome file, integrate
				String integChrFileName = integBaseName + spacer + "integ_" + chrName + ".bedGraph";
				Callable<File> worker = new FileIntegration(matchChrFiles, integChrFileName);
				Future<File> submit = executor.submit(worker);
				futureList.add(submit);
			}
			//else skip integration
		}

		for (Future<File> future : futureList) {
			try {
				File f = future.get();
				System.out.println(f.getPath() + " integration completed");
				outFileList.add(f);
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
		}

		executor.shutdown();
		executor.awaitTermination(12, TimeUnit.HOURS);

		return;
	}

	/**
	 * Given a file f, parse the file name for the chromosome name
	 * @param f
	 * @return chromosome name - (i.e. chr17)
	 */
	private static String parseBGChrom(File f) {
		String filename = f.getName();
		return filename.substring(filename.lastIndexOf("_chr") + 1, filename.lastIndexOf(".bedGraph"));
	}

}
