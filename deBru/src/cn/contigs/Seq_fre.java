package cn.contigs;

import java.util.Objects;

public class Seq_fre implements Comparable<Seq_fre>{
    public long seq;
    public Long fre;

    public Seq_fre(){}
    public Seq_fre(long seq, long fre) {
        this.seq = seq;
        this.fre = fre;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Seq_fre seq_fre = (Seq_fre) o;
        return seq == seq_fre.seq;
    }

    @Override
    public int hashCode() {

        return Objects.hash(seq);
    }

    @Override
    public int compareTo(Seq_fre o) {
        return this.fre.compareTo(o.fre);
    }
}
