
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;
import java.util.concurrent.Callable;

/**
 * Integrates the given bedGraph files
 * Output goes to first filename
 */
public class FileIntegration implements Callable<File> {
	
	/** Used when calling from main */
	private static ArrayList<String> filenames;
	/** Used when calling from GMIM Pipeline */
	private static ArrayList<File> files;
	/** Integration output filename */
	private static String outfile;
	
	public FileIntegration(ArrayList<File> fs, String oF) {
		files = fs;
		outfile = oF;
	}
	
	/**
	 * Main method
	 *@param args - bedGraph files to integrate
	 */
	public static void main(String[] args) throws Exception{
		filenames = new ArrayList<String>();
		
		if (args.length < 3) {
			System.err.println("Must include files to integrate");
			System.exit(1);
		}
		
		outfile = args[0];
		
		for (int i = 1; i < args.length; i++) {
			String s = args[i];
			if (!s.endsWith(".bedGraph")) {
				System.err.println("Incorrect filetypes - must be a .bedGraph: " + s);
				System.exit(1);
			}
			
			filenames.add(s);
		}
		
		if (filenames.isEmpty()) {
			System.err.println("No files to integrate");
			System.exit(1);
		}
		
		startIntegMain();
	}
	
	@Override
	public File call() throws Exception {
		ArrayList<Scanner> scannerList = new ArrayList<Scanner>();

		for (File f: files) {
			try {
				Scanner sc = new Scanner(f);
				scannerList.add(sc);
				sc.nextLine(); //ignore header
			} catch (FileNotFoundException e) {
				System.err.println("File not found: " + f.getPath());
				System.exit(1);
			}
		}
		
		File integFile = integration(scannerList);
		return integFile;
		
	}

	public static void startIntegMain() throws FileNotFoundException {
		ArrayList<Scanner> scannerList = new ArrayList<Scanner>();

		for (String filename : filenames) {
			try {
				Scanner sc = new Scanner(new File(filename));
				scannerList.add(sc);
				sc.nextLine(); //ignore header
			} catch (FileNotFoundException e) {
				System.err.println("File not found: " + filename);
				System.exit(1);
			}
		}
		
		integration(scannerList);
	}
	
	/**
	 * Integrate the output bedGraph files
	 * @throws FileNotFoundException 
	 * @throws IOException 
	 */
	public static File integration(ArrayList<Scanner> scannerList) throws FileNotFoundException {
		
		File integFile = new File(outfile);
		PrintWriter pw = null;
		try {
			
			pw = new PrintWriter(integFile);
		} catch (FileNotFoundException e) {
			throw new FileNotFoundException("File not found: " + outfile);
		}
		
		boolean endReached = false;
		
		pw.println("track type=bedGraph name=\"" + outfile.substring(0, outfile.lastIndexOf(".bedGraph")) + "\" " 
				+ "visibility=full autoScale=Off alwaysZero=On maxHeightPixels=128:30:11 viewLimits=0:1 color=255,30,30");
		
		while(scannerList.get(0).hasNextLine() && endReached == false) {
			String chr = "";
			int start = -1;
			int end = -1;
			
			double product = 1;
			String[] raComp = null;
			
			for (int i = 0; i < scannerList.size() && !endReached; i++) {
				Scanner sc = scannerList.get(i);
				
				if (!sc.hasNextLine()) {
					endReached = true;
				}
				else { //process file
					
					String line = sc.nextLine();
					String[] ra = line.split("\\s");
						
					//check that the line chr/start/end match
					if (raComp == null) {
						raComp = Arrays.copyOf(ra, ra.length);
					}
					else if (!lineMatchingRegion(raComp, ra)) {
						System.err.println("File regions do not match: " + line);
						System.exit(1);
					}
					
					chr = ra[0];
					start = Integer.parseInt(ra[1]);
					end = Integer.parseInt(ra[2]);			
					double cMBF = Double.parseDouble(ra[3]);
					product *= cMBF;
				}
			}
			pw.println(chr + "\t" + start + "\t" + end + "\t" + product);
		}

		for (Scanner sc: scannerList) {
			sc.close();
		}
		pw.close();
		
		return integFile;
	}

	private static boolean lineMatchingRegion(String[] raComp, String[] ra) {
		for (int i = 0; i < 3; i++) { //check chr, start, end
			if (!raComp[i].equals(ra[i])) {
				return false;
			}
		}
		return true;
	}
	

}
