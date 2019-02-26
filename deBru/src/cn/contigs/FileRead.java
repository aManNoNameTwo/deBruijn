package cn.contigs;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;
import java.io.File;
import java.util.List;
import java.util.Map;

public class FileRead {
    public static Map<String, Object> BufferedRandomAccessFileReadLine(File file, String encoding, long pos, int num) {
        Map<String, Object> res = Maps.newHashMap();
        List<String> pins = Lists.newArrayList();
        res.put("pins", pins);
        BufferedRandomAccessFile reader = null;
        try {
            reader = new BufferedRandomAccessFile(file, "r");
            reader.seek(pos);
            for (int i = 0; i < num; i++) {
                String pin = reader.readLine();
                if (org.apache.commons.lang3.StringUtils.isBlank(pin)) {
                    break;
                }
                pins.add(new String(pin.getBytes("8859_1"), encoding));
            }
            res.put("pos", reader.getFilePointer());
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            IOUtils.closeQuietly(reader);
        }
        return res;
    }
}
