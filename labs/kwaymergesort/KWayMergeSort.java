package kwaymergesort;

import timing.Ticker;

import java.util.ArrayList;
import java.util.LinkedList;
public class KWayMergeSort {
	
	/**
	 * 
	 * @param K some positive power of 2.
	 * @param input an array of unsorted integers.  Its size is either 1, or some other power of 2 that is at least K
	 * @param ticker call .tick() on this to account for the work you do
	 * @return
	 */
	public static Integer[] kwaymergesort(int K, Integer[] input, Ticker ticker) {
		int n = input.length;
		
		//
		// FIXME
		// Following just copies the input as the answer
		//
		// You must replace the loop below with code that performs
		// a K-way merge sort, placing the result in ans
		//
		// The web page for this assignment provides more detail.
		//
		// Use the ticker as you normally would, to account for
		// the operations taken to perform the K-way merge sort.
		//
		//if size = 0 or 1, return input
		Integer[] ans = new Integer[n];
		if(n <= 1){
			ticker.tick();
			return input;
		}
		//else recursively call merge sort on list  (left:K, right n-K)
		else{
			Integer[][] L = new Integer[K][n/K];
			for(int i = 0; i < n; i++){
				L[i/(n/K)][i%(n/K)] = input[i];		
				ticker.tick();
			}
			for(int i = 0; i < K; i++) {
				L[i] = kwaymergesort(K, L[i], ticker);
				ticker.tick();
			}
			ans = merge(L, ticker);
			ticker.tick();
		}
		return ans;
	}
	public static Integer[] merge(Integer[][] pieces, Ticker ticker){
        int k = pieces.length;
        int size = pieces[0].length;
        if (k == 1) {
        	ticker.tick();
        	return pieces[0];        
        }
        else {
        int[] track = new int[k];
        Integer[] merged = new Integer[k * size];
        int counter = 0;
        ticker.tick();
        while (counter < k * size){
            int min = Integer.MAX_VALUE;
            int checking = 0;
            ticker.tick(); 
            for (int i = 0; i < k; i++){
                if (track[i] < size && pieces[i][track[i]] < min){
                        min = pieces[i][track[i]];
                        checking = i;
                        ticker.tick();
                    }
                }                
            track[checking] = track[checking]+1;            
            merged[counter] = min;
            counter = counter +1;		
            ticker.tick();
    	}
        return merged;
        }

        }
    }
