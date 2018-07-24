package wGSepPar_CLI;
import java.util.ArrayList;

public class IntWindow_WGSep {
	private static int intSize;
	private static double medianMult;
	
	/** window size in terms of lines **/
	private int windowSize;
	
	private String chromNum;
	private int absStart;
	private int start;
	private int pos;
	private int end;
	private int relIndex;
	private boolean endOfChrom;

	protected ArrayList<IntStats_WGSep> winRA;

	/**
	 * @param wS - window size in terms of lines
	 * @param iS - interval size
	 */
	public IntWindow_WGSep(String cN, int wS, int iS, double mM, int p) {
		chromNum = cN;
		intSize = iS;
		windowSize = wS / intSize;
		winRA = new ArrayList<IntStats_WGSep>(windowSize);
		
		medianMult = mM;
		endOfChrom = false;
		pos = p;
		
		absStart = p;
		start = p;
		relIndex = 0;
		//calcStartEnd();
	}

	private void calcStartEnd() {
		if (!endOfChrom) { //don't change window start/end if reached end of file
			if (windowSize % 2 != 0) { //if window size odd
				start = pos - intSize * (windowSize/ 2);
				end = pos + intSize * ((windowSize + 1) / 2);
			}
			else {
				start = (pos - (int) (Math.floor((windowSize - 1)/2) * intSize));
				end = (pos + (int) (Math.ceil((windowSize + 1)/2) * intSize));
			}
		}
		
		if (start < absStart) {
			end = (end - start) + absStart;
			start = absStart;
		}
		
		relIndex = (pos - start) / intSize;
	}

	public void setEndOfChrom() {
		endOfChrom = true;
	}
	
	//Used for when window is not filled but reached end of file - assumingly 0 - [last line's start]
	public int setSmallerWindowSize() {
		windowSize = winRA.size();
		end = winRA.get(winRA.size() - 1).getEnd();
		return windowSize;
	}

	public void setChromNum(String cN) {
		chromNum = cN;
	}
	
	public String getChromNum() {
		return chromNum;
	}
	
	public int getIndexStart() {
		return winRA.get(relIndex).getStart();
	}

	public int getIndexEnd() {
		return winRA.get(relIndex).getEnd();
	}

	public double getIndexReadCount() {
		return winRA.get(relIndex).getReadCount();
	}

	public double getIndexcMBF() {
		return winRA.get(relIndex).getcMBF();
	}
	
	public int getLastStartIndex() {
		return winRA.get(winRA.size() - 1).getStart();
	}

	/**
	 * Calculates the number of elements left to fill in window
	 * @return
	 */
	public int toFill() {
		return windowSize - winRA.size();		
	}

	public boolean full() {
		return winRA.size() == windowSize;
	}

	public boolean insert(IntStats_WGSep ps) {
//		if (!ps.getChromNum().equals(chromNum)) {
//			throw new IllegalArgumentException("Cannot add a position from a different chromosome");
//		}
//		else 
		if (winRA.size() >= windowSize) { //full
			throw new IndexOutOfBoundsException();
		}
		return winRA.add(ps);
	}

	private double calcMean() { //used in place of Z
		double sum = 0;
		int count = 0;
		for (int i = 0; i < windowSize; i++) {
			sum += winRA.get(i).getReadCount();
			count++;
		}
		return sum / count;
	}

	private double calcMedian() {
		@SuppressWarnings("unchecked")
		ArrayList<IntStats_WGSep> sortedRA = (ArrayList<IntStats_WGSep>) winRA.clone();
		sortedRA.sort(new ISWGSepReadCountSorter());
		int middle = 0;
		if (windowSize % 2 == 0) { //if even, average middle two entries
			middle = windowSize/2 - 1;
			return (sortedRA.get(middle).getReadCount() + sortedRA.get(middle + 1).getReadCount()) / 2;
		}
		else {
			middle = windowSize/2;
			return sortedRA.get(middle).getReadCount(); //else odd, return middle entry
		}
		
	}

	public double calccMBF() throws Exception{ //currently uses mean
		double noiseEst = medianMult * calcMedian();
		if (noiseEst == 0) {
			throw new Exception("Noise level estimated to be 0: \n"
					+ chromNum + "\t" + start + "\t" + end + "\t" + pos);
		}
		
		double Z = getIndexReadCount() / noiseEst;
		double cmbf = 1 - Math.exp(-1 * (Z*Z) / 2);
		winRA.get(relIndex).setcMBF(cmbf);
		return cmbf;
	}

	public void incrCenter() {
		pos += intSize;
		calcStartEnd();
		if (!endOfChrom && winRA.get(0).getStart() < start) { //remove first element in RA if moving out of range and haven't reached end of file
			winRA.remove(0);
		}
	}



	@Override
	public String toString() {
		return "windowSize:" + windowSize + ", pos:" + pos + ", start:" + start + ", end:" + end + ", relIndex" + relIndex + "\nwinRA:" + winRA.toString();
	}

}
