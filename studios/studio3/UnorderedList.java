package studio3;

import java.util.LinkedList;
import java.util.List;

public class UnorderedList<T extends Comparable<T>> implements PriorityQueue<T> {

	public List<T> list;
	
	public UnorderedList() {
		list = new LinkedList<T>();
	}
	
	@Override
	public boolean isEmpty() {
		if (list.size() == 0) {
		return false;
		}
		else {
		return true;
		}	
	}

	@Override
	public void insert(T thing) {
		int header;
		int e;
		
		
		addBefore(e,header);
	}

	@Override
	public T extractMin() {
		//
		// FIXME
		//
		return null;
	}

}
