package wGSepPar_CLI;
import java.util.Comparator;

/**
 * Object that holds the read count and cMBF of a position in the chromosome
 * Can get/set readCount/cMBF
 */
class IntStats_WGSep {
	
	private static double minRC = 0.1;
	
	private String _chromNum;
	private int _start;
	private int _end;
	private int _width;
	private double _readCount;
	private double _cMBF;
	
	
	public IntStats_WGSep(String cN, int s, int e, int w, int rc) {
		_chromNum = cN;
		_start = s;
		_end = e;
		_width = w;

		setReadCount(rc);
		
		_cMBF = -1;
	}
	
	public static void setMinRC(double mRC) {
		minRC = mRC;
	}
	
	public String getChromNum() {
		return _chromNum;
	}
	
	public int getStart() {
		return _start;
	}
	
	public int getWidth() {
		return _width;
	}
	
	public int getEnd() {
		return _end;
	}
	
	public double getReadCount() {
		return _readCount;
	}
	
	public double getcMBF() {
		return _cMBF;
	}
	
	public void setReadCount(int rc) {
		if (rc > 0) {
			_readCount = rc;
		}
		else if (rc == 0){
			_readCount = minRC;
		}
		else {
			System.err.println("readCount cannot be negative" + rc);
			System.exit(1);
		}
	}
	
	
	public void setcMBF(double cmbf) throws Exception {
		if (cmbf >= 0 && cmbf <= 1) {
			_cMBF = cmbf;
		}
		else{
			throw new Exception("cMBF cannot be outside the range of 0 to 1: " + cmbf + "\n"
					+ this._chromNum + "\t" + this._start + "\t" + this._end + "\t" + this._readCount);
		}
	}
	
	@Override
	public String toString() {
		return "[start:" + _start + ", end:" + _end + ", readCount:" + _readCount + ", cMBF:" + _cMBF + "]";
	}
	
}

class ISWGSepReadCountSorter implements Comparator<IntStats_WGSep> {

	@Override
	public int compare(IntStats_WGSep ps1, IntStats_WGSep ps2) {
		return (int) (10*(ps1.getReadCount() - ps2.getReadCount()));
	}
	
}