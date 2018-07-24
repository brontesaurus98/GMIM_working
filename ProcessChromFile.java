package wGSepPar_CLI;

import java.io.File;
import java.io.PrintWriter;
import java.util.Scanner;
import java.util.concurrent.Callable;

public class ProcessChromFile implements Callable<File> {
	/** Window size for calculating cMBF (in bp) */
	private static Integer windowbpSize;
	/** Read resolution / width of bin */
	private static int intervalSize;
	/** Median multiple */
	private static double medianMult;
	private static String outBaseName = null;

	private String dirName;
	private File chromFile;

	public ProcessChromFile(String dN, File cF) {
		dirName = dN;
		chromFile = cF;
	}
	
	public static void setWinSize(int wS) {
		windowbpSize = wS;
	}
	
	public static void setIntSize(int iS) {
		intervalSize = iS;
	}
	
	public static void setMedMult(double mM) {
		medianMult = mM;
	}
	
	public static void setOutBaseName(String oBN) {
		outBaseName = oBN;
	}

	@Override
	public File call() throws Exception {
		String outDirName = dirName.substring(0, dirName.lastIndexOf("_ChromFiles")) + "_out";
		String spacer = "/";
		File outDir = new File(outDirName);
		boolean outFolderMade = outDir.mkdir(); 
		if (!outDir.exists() && !outFolderMade) {
			spacer = "_";
		}
		
		String outFileName = "";
		String chromFileName = chromFile.getName().substring(0, chromFile.getName().length() - 4);
		
		if (outBaseName == null) {
			outFileName = "out_" + chromFileName;
		}
		else {
			outFileName = outBaseName + "_" + chromFileName;
		}
		
		String outFilePath = outDirName + spacer + outFileName;
		File outFile = new File(outFilePath + ".bedGraph");
		PrintWriter pw = new PrintWriter(outFile);
		pw.println("track type=bedGraph name=\"" + outFilePath + "\"" + " description=\"" + outFilePath + "\" "
				+ "visibility=full autoScale=Off alwaysZero=On maxHeightPixels=128:30:11 viewLimits=0:1"); //header

		Scanner sc = new Scanner(chromFile);
		if (!sc.hasNextLine()) {
			sc.close();
			pw.close();
			return outFile;
		}

		String line = sc.nextLine();
		IntStats_WGSep startIS = parseLine(line);
		
		if (startIS == null) {
			System.err.println("Unable to parse line: " + line);
			System.exit(1);
		}

		int index = startIS.getStart();
		intervalSize = startIS.getEnd() - startIS.getStart();
		String chromNum = startIS.getChromNum();
		
		if (windowbpSize % intervalSize != 0) {
			sc.close();
			pw.close();
			System.err.println("Window size must be a multiple of the interval size");
			System.exit(1);
		}

		IntWindow_WGSep window = new IntWindow_WGSep(chromNum, windowbpSize, intervalSize, medianMult, index);
		window.insert(startIS);

		//parse through file for position index and its read count
		while(sc.hasNextLine()) {
			if (window.toFill() > 0) {
				line = sc.nextLine();
				IntStats_WGSep ps = parseLine(line);
				if (ps == null) continue;
				window.insert(ps);
			}

			if (window.full()) { //filled window, calculate stats for current position and increment (middle indices)
				printIndexStats(pw, window);

				if (!sc.hasNextLine()) {
					window.setEndOfChrom();
				}

				window.incrCenter();
				index += intervalSize;
			}
		}

		if (!sc.hasNextLine()) { //reached end of file, finish computing for last indices
			window.setEndOfChrom();

			if (!window.full()) { //reached end of file, but window not filled
				windowbpSize = window.setSmallerWindowSize() * intervalSize;
			}

			int lastIndex =  window.getLastStartIndex();

			for (index += 0; index <= lastIndex; index += intervalSize) {
				printIndexStats(pw, window);
				window.incrCenter();
			}

		}
		
		sc.close();
		pw.close();
		return outFile;
	}

	public IntStats_WGSep parseLine(String line) {
		String[] ra = line.split("\\s+");
		if (line.isEmpty() || ra.length == 0) return null;

		try {
			String chromNum;
			int intStart, intEnd, intReadCount;
			
			chromNum = ra[0];
			intStart = Integer.parseInt(ra[1]);
			intEnd = Integer.parseInt(ra[2]);
			intReadCount = Integer.parseInt(ra[3]);
			return new IntStats_WGSep(chromNum, intStart, intEnd, intervalSize, intReadCount);
		} catch (NumberFormatException e) {
			System.err.println("File has improper values: " + line);
			System.exit(1);
		} catch (IndexOutOfBoundsException e) {
			System.err.println("Improper file formatting/values: " + line);
			System.exit(1);
		}

		return null;
	}

	public void printIndexStats(PrintWriter pw, IntWindow_WGSep window) {
		try {
			String cN = window.getChromNum();
			int iS = window.getIndexStart();
			int iE = window.getIndexEnd();
			double icMBF = window.calccMBF();

			pw.printf("%s\t%d\t%d\t%.5f\n", cN, iS, iE, icMBF);
		} catch (Exception e) {
			System.err.println(e.getMessage());
			System.exit(1);
		}
	}
}
