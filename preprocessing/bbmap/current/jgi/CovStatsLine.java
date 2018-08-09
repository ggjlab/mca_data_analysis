package jgi;

import java.util.Arrays;

import align2.Tools;

/**
 * @author Brian Bushnell
 * @date Oct 10, 2014
 *
 */
public class CovStatsLine {
	
	public CovStatsLine(String s){
		this(s.split("\t"));
	}

	/**
	 * ID	Avg_fold	Length	Ref_GC	Covered_percent	Covered_bases	Plus_reads	Minus_reads	(optional.... Read_GC)
	 * @param split
	 */
	public CovStatsLine(String[] split) {
		assert(split.length>=8) : Arrays.toString(split);
		assert(!split[0].startsWith("#")) : Arrays.toString(split);
		id=split[0];
		avgFold=Double.parseDouble(split[1]);
		length=Integer.parseInt(split[2]);
		refGC=Double.parseDouble(split[3]);
//		coveredPercent=Double.parseDouble(split[4]);
		coveredBases=Integer.parseInt(split[5]);
		plusReads=Long.parseLong(split[6]);
		minusReads=Long.parseLong(split[7]);
		if(split.length==11){
			median=Integer.parseInt(split[8]);
			underMin=Integer.parseInt(split[9]);
			readGC=Double.parseDouble(split[10]);
		}else if(split.length==10){
			median=Integer.parseInt(split[8]);
			if(CoveragePileup.USE_BITSETS && CoveragePileup.USE_WINDOW){
				underMin=Integer.parseInt(split[9]);
			}else{
				readGC=Double.parseDouble(split[9]);
			}
		}else if(split.length==9){
			readGC=Double.parseDouble(split[8]);
		}else if(split.length<9){
			//do nothing
		}
	}
	
	public final double coveredPercent(){
		return (100.0*coveredBases)/Tools.max(1, length);
	}
	
	public final long reads(){return plusReads+minusReads;}
	
	/**
	 * @param csl
	 */
	public void add(CovStatsLine csl) {
		double invlen2=1d/Tools.max(1, length+csl.length);
		avgFold=((avgFold*length)+(csl.avgFold*csl.length))*invlen2;
		refGC=((refGC*length)+(csl.refGC*csl.length))*invlen2;
		readGC=((readGC*reads())+(csl.readGC*csl.reads()))*1.0/(Tools.max(1, reads()+csl.reads()));
		
		length+=csl.length;
		coveredBases+=csl.coveredBases;
		plusReads+=csl.plusReads;
		minusReads+=csl.minusReads;
		median=median+csl.median;
		underMin=underMin+csl.underMin;
	}
	
	public String toString(){
		return String.format("%s\t%.4f\t%d \t%.4f\t%.4f\t%d \t%d\t%d\t%d\t%d\t%.4f", id, avgFold, length,
				refGC, coveredPercent(), coveredBases, plusReads, minusReads, median, underMin, readGC);
	}
	
	public static final String header1="#ID\tAvg_fold\tLength\tRef_GC\tCovered_percent\tCovered_bases\tPlus_reads\tMinus_reads\tMedian_fold\tUnder_min\tRead_GC";
	public static final String header2="#ID\tAvg_fold\tLength\tRef_GC\tCovered_percent\tCovered_bases\tPlus_reads\tMinus_reads\tRead_GC";
	
	public String id;
	public int length;
	public int coveredBases;
	public long plusReads;
	public long minusReads;
	public double avgFold;
	public double refGC;
	public int median;
	public int underMin;
	public double readGC;
}
