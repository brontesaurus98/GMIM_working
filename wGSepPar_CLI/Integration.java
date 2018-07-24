import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

/**
 * Integrates the given bedGraph files
 * Output goes to first filename
 */
public class Integration {
	
	private static ArrayList<String> filenames;
	private static String outfile;
	
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
		
		integration();
	}
	
	/**
	 * Integrate the output bedGraph files
	 * @throws IOException 
	 */
	public static void integration() {

		ArrayList<Scanner> scannerList = new ArrayList<Scanner>();
		PrintWriter pw = null;
		
		try {
			pw = new PrintWriter(new File(outfile));
		} catch (FileNotFoundException e) {
			
		}
		
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
		
		boolean endReached = false;
		
		pw.println("track type=bedGraph name=\"" + outfile.substring(0, outfile.lastIndexOf(".bedGraph")) + "\"");
		
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
						System.err.println("File regions do not match" + line);
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
