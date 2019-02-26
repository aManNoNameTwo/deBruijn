package cn.contigs;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

public class ReadID_Pos {
    public String readID;
    public List<Integer> readPos = new ArrayList<>();
    public boolean deleFlag;

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        ReadID_Pos that = (ReadID_Pos) o;
        return Objects.equals(readID, that.readID);
    }

    @Override
    public int hashCode() {

        return Objects.hash(readID);
    }
}
