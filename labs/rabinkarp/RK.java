package rabinkarp;

public class RK {
	
	//
	// Be sure to look at the write up for this assignment
	//  so that you get full credit by satisfying all
	//  of its requirements
	//
	
	/**
	 * Rabin-Karp string matching for a window of the specified size
	 * @param m size of the window
	 */
	public int window_size;
	public char[] window;
	public int hash;
	public int hash_m;
	public int l = 0;
	public int mod(int m) {
		   int hash_m = 1;
		   for (int i = 1; i <= m;) {
				hash_m = (31*hash_m)%511;
				i = i+1;
			}
		   return hash_m;
	}
	 
	public RK(int m) {
		this.window_size = m;
		this.window = new char[this.window_size];
		this.hash = 0;
	 }
	/**
	 * Compute the rolling hash for the previous m-1 characters with d appended.
	 * @param d the next character in the target string
	 * @return
	 *
	 */


	public int nextCh(char d) {
		int m = window.length;
		int c = window[l];
		window[l]=d;
		l = l+1;
		if(window[window_size-1]==0){
			hash = (hash * 31 + window[l-1]  ) % 511;
		}else{
			hash = (((hash * 31) - (mod(m)*c) + d) % 511);
			l = l%window_size;
		}
		if(hash<0){
			hash = hash + 511;
		}
		return hash;
	}
	
	


}