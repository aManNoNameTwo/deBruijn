package cn.contigs;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

public class K_mers {
    private long base_seq;
    //private K_mers nextK_mers;
    private long numOfRead;
    private List<ReadID_Pos> readID_Pos;
    private long frequence;
    public int num_visit =0;

    public long getNumOfRead() {
        return numOfRead;
    }

    public void setNumOfRead(long numOfRead) {
        this.numOfRead = numOfRead;
    }
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        K_mers k_mers = (K_mers) o;
        return Objects.equals(base_seq, k_mers.base_seq);
    }

    @Override
    public int hashCode() {

        return Objects.hash(base_seq);
    }


    public long getFrequence() {
        return frequence;
    }

    public void setFrequence(long frequence) {
        this.frequence = frequence;
    }

    public K_mers(){}

    public K_mers(long base_seq, List<ReadID_Pos> readID_Pos, long frequence,long numOfRead) {
        this.base_seq = base_seq;
        //this.nextK_mers = nextK_mers;
        this.readID_Pos = readID_Pos;
        this.frequence = frequence;
        this.numOfRead = numOfRead;
    }
    public long getbase_seq() {
        return base_seq;
    }

    public void setbase_seq(long base_seq) {
        this.base_seq = base_seq;
    }

//    public K_mers getNextK_mers() {
//        return nextK_mers;
//    }
//
//    public void setNextK_mers(K_mers nextK_mers) {
//        this.nextK_mers = nextK_mers;
//    }

    public List<ReadID_Pos> getReadID_Pos() {
        return readID_Pos;
    }

    public void setReadID_Pos(List<ReadID_Pos> readID_Pos) {
        this.readID_Pos = readID_Pos;
    }
}
