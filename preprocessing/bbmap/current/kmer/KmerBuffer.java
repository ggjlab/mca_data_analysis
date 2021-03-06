package kmer;

import stream.ByteBuilder;
import align2.IntList;
import align2.LongList;

/**
 * @author Brian Bushnell
 * @date Jul 30, 2015
 *
 */
public class KmerBuffer {
	
	public KmerBuffer(int buflen, int k_, boolean initValues){
		k=k_;
		kmers=new LongList(buflen);
		values=(initValues ? new IntList(buflen) : null);
	}
	
	public int add(long kmer){
		assert(values==null);
		kmers.add(kmer);
		return kmers.size;
	}
	
	public void add(long kmer, int value){
		kmers.add(kmer);
		values.add(value);
		assert(values.size==kmers.size);
	}
	
	public void clear(){
		kmers.clear();
		if(values!=null){values.clear();}
	}
	
	//Returns raw size of kmers array, rather than actual number of kmers
	final int size(){return kmers.size;}
	
	public String toString(){
		ByteBuilder bb=new ByteBuilder();
		for(int i=0; i<kmers.size; i++){
			if(i>0){bb.append(',');}
			bb.appendKmer(kmers.get(i), k);
		}
		return bb.toString();
	}
	
	private final int k;
	final LongList kmers;
	final IntList values;
	
}
