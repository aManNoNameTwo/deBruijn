package cn.contigs;
import scala.Int;

import java.io.*;
import java.util.*;
import java.util.function.IntToDoubleFunction;
import java.util.regex.Pattern;

public class ContigsCreation implements Constants{
    public static int k = 11  ;
    public static int deep =60;
    public static float cur_rate =0.6f;
    public long addSum = sumOf_2(2*k);
    public Seq_fre seq_fre = new Seq_fre(0,0);
    public Map<Long,K_mers> kmersHash = new HashMap();
    public DeBruijnGraph debru = new DeBruijnGraph();

    public static void main(String[] args){
        //String filepath = "D:\\SRR835425.fastq";
        String filepath = "D:\\chr1_197448179_197448200_18_31.fa";
        ContigsCreation conc = new ContigsCreation();
        conc.DeBruijn(filepath); //start

    }
    public String DeBruijn(String filePath){

        //1.get reads and kmers
        try{
            kmersHash = readFile_getKmer(filePath);
        }catch (Exception e){
            e.printStackTrace();
            return null;
        }
        K_mers firstKmer = createInitilKmer(kmersHash);
        debru.addVex(firstKmer);
        createDebru(firstKmer);
        if(kmersHash.containsKey(bit_Reverse(firstKmer.getbase_seq()))){
            createDebru(kmersHash.get(bit_Reverse(firstKmer.getbase_seq())));
        }
        //debru.showGraph("preGraph");
        merge_unipath();

        debru.printf_reEdge();
        debru.showGraph("afterGraph");
        return null;
    }
    public long bit_Reverse(long v) {   //reverse the number
        v = ((v>>> 1) & 0x5555555555555555L) | ((v << 1) & 0xaaaaaaaaaaaaaaaaL);
        v = ((v >>> 2) & 0x3333333333333333L) | ((v << 2) & 0xccccccccccccccccL);
        v = ((v >>> 4) & 0x0f0f0f0f0f0f0f0fL) | ((v << 4) & 0xf0f0f0f0f0f0f0f0L);
        v = ((v >>> 8) & 0x00ff00ff00ff00ffL) | ((v << 8) & 0xff00ff00ff00ff00L);
        v = ((v >>> 16) & 0x0000ffff0000ffffL) | ((v << 16) & 0xffff0000ffff0000L);
        v = ((v >>> 32) & 0x00000000ffffffffL) | ((v << 32) & 0xffffffff00000000L);
        v = v >>> (64-2*k);
        return v;
    }
    public boolean createDebru(K_mers kmer){
        int cur_deep = 0;
        boolean isSuccessEdge = true;
        K_mers pre_kmer = new K_mers();
        K_mers next_kmer = new K_mers();
        while(kmer != null && cur_deep <=deep){

            next_kmer = greedy_nextKmer(kmer.getbase_seq(),kmer.num_visit%4);
            if(next_kmer != null){
                debru.addVex(next_kmer);
                isSuccessEdge = debru.addEdge(kmer,next_kmer);
                if(!isSuccessEdge){
                    kmer.num_visit ++;
                    continue;
                }else{
                    cur_deep ++;
                    pre_kmer = kmer;
                    kmer = next_kmer;

                }

            }else {
                long position = debru.getPosition(kmer,debru.vertexList);
                if(debru.vertexList.get((int)position).firstIn == null)    break;
                kmer = debru.vertexList.get((int)debru.vertexList.get((int)position).firstIn.tailvex).kmer;
                kmer.num_visit ++;
            }
            if(kmer.num_visit >= 4){
                long pos = 0L;
                boolean is_end = false;
                while(kmer.num_visit >= 4){  // loop for find the kmer which have small num_visit
                    pos = debru.getPosition(kmer,debru.vertexList);
                    if(debru.vertexList.get((int)pos).firstIn == null  ||
                            debru.vertexList.get((int)pos).firstIn.tailvex == pos){
                        is_end = true;
                        break;
                    }
                    kmer = debru.vertexList.get((int)debru.vertexList.get((int)pos).firstIn.tailvex).kmer;
                    kmer.num_visit ++;
                }
                if(is_end) break;
            }

        }
        return true;
    }
    public boolean createDebru(K_mers kmer,int ori,int cur_deep){
        if(kmer == null) return false;
        if(cur_deep == deep) return false;

        if(ori == 0){
            K_mers next_kmer = greedy_nextKmer(kmer.getbase_seq(),1);
            if(next_kmer == null) return false;
            debru.addVex(next_kmer);
            if(!debru.addEdge(kmer,next_kmer)){
                return false;
            }
            //kmersHash.remove(next_kmer.getbase_seq());
            createDebru(next_kmer,0,cur_deep+1);
            createDebru(kmer,0,cur_deep);
        }else if(ori ==1){
            long reverse_Kmer =bit_Reverse(kmer.getbase_seq());
            if(kmersHash.containsKey(reverse_Kmer)){
                K_mers next_kmer = greedy_nextKmer(reverse_Kmer,1);
                if(next_kmer == null) return false;
                debru.addVex(next_kmer);
                if(!debru.addEdge(kmer,next_kmer)){
                    return false;
                }
                //kmersHash.remove(next_kmer.getbase_seq());
                createDebru(next_kmer,0,cur_deep+1);
                createDebru(kmer,0,cur_deep);
            }
        }
        return true;
    }
    public K_mers greedy_nextKmer(long base_Seq,int kth){
        //System.out.println(num_visit);
        long seq = base_Seq;
        base_Seq = (base_Seq << 2) & (addSum);
        long[] fre_seq = new long[4];
        for(int i=0; i< 4;i++){
            long addre = base_Seq|i; //get the base_seq next seq
            if(kmersHash.containsKey(addre) && seq != addre){
                K_mers cur_kmer = kmersHash.get(addre);
                long cur_fre = cur_kmer.getFrequence();
                if(cur_fre >= MIN_FREQUENCE){
                    fre_seq[i] = cur_fre;
                }else{
                    fre_seq[i] = -1L;
                }
            }else{
                fre_seq[i] = -1L;
            }
        }
        long k_max = -1L;
        long[] fre_seqTem = fre_seq.clone();
        k_max = findKth(fre_seqTem,kth+1);
        int indexmax = -1;
        if(k_max ==-1){
            return null;
        }else{
            for(int i=0;i<fre_seq.length;i++){
                if(fre_seq[i] == k_max){
                    indexmax = i;
                    break;
                }
            }
            return kmersHash.get(base_Seq |indexmax);
        }
    }
    public long findKth(long[] array,int k){
        Arrays.sort(array);
        //System.out.printf("length:%d -------k:%d\n",array.length,k);
        return array[array.length-k];
    }
    public  long findKth(long[] a, int n, int K) {
        return findKth(a, 0, n - 1, K);
    }
    public  long findKth(long[] a, int start, int end, int k) {
        int pivot = partation(a, start, end);
        if (k == pivot - start + 1){
            return a[pivot];
        } else if (k < pivot - start + 1) {
            return findKth(a, start, pivot - 1, k);
        }else {
            return findKth(a,pivot+1,end,k-(pivot - start +1));
        }
    }
    public  int partation(long[] a, int low, int high) {
        long key = a[low];
        while (low < high) {
            while (low < high && a[high] <= key)
                high--;
            a[low] = a[high];
            while (low < high && a[low] >= key)
                low++;
            a[high] = a[low];
        }
        a[low] = key;
        return low;
    }
    public Map<Long,K_mers> readFile_getKmer (String filePath)throws Exception {
        String filetype = filePath.split("\\.")[1];
        if(Pattern.matches("^f(ast)?a$",filetype)){
            filetype = "fasta";
            return readFa(filePath);
        }else if(Pattern.matches("^f(ast)?q$",filetype)){
            filetype = "fastq";
            return readFa(filePath);
        }else {
            throw new Exception("only can read fastq or fasta file");
        }
    }
    public Map<Long, K_mers> readFq(String filePath){
        final String ENCODING = "UTF-8";
        final int NUM = 4000;
        Map<Long,K_mers> kmersHash = new HashMap();  //k-mer s库
        File file = new File(filePath);
        long pos = 0L;
        String M ="A";
        int nu =0;
        while (true) {
            Map<String, Object> res = FileRead.BufferedRandomAccessFileReadLine(file, ENCODING, pos, NUM);
            // 如果返回结果为空结束循环
            if (org.apache.commons.collections4.MapUtils.isEmpty(res)) {
                break;
            }
            List<String> pins = (List<String>) res.get("pins");
            if (org.apache.commons.collections4.CollectionUtils.isNotEmpty(pins)) {
                //generate kmers
                for(int i=1; i<pins.size();i+=4){
                    //System.out.println(pins.get(i).length());
                    if(pins.get(i).contains("N")){
                        continue;
                    }
                    int count = (pins.get(i).length()-pins.get(i).replace(M, "").length())/M.length();
                    if(count/pins.get(i).length() > 0.8){
                        continue;
                    }
                    String qual = pins.get(i+2);
                    float sum_cur_rate = 0f;
                    for(int j=0; j<qual.length();j++){
                        sum_cur_rate = sum_cur_rate + (float)Math.pow(10.0,-((double)((qual.charAt(j)-33)/10.0)));  //calculate
                    }
                    if((1.0-sum_cur_rate/(float)qual.length()) >= cur_rate){
                        //System.out.println(nu++);
                        String read_seq = pins.get(i);  //get the right read
                        long k_mer = 0;
                        if(k < read_seq.length()){
                            for(int L=0; L<read_seq.length()-k+1;L++){
                                if(L==0){
                                    k_mer = kmersAdders(read_seq.substring(0,k));
                                    //System.out.println(Long.toBinaryString(k_mer));
                                }else{
                                    k_mer =((k_mer <<2) | pre_base(read_seq.charAt(L+k-1)))&(addSum);
                                    //System.out.println(Long.toBinaryString(k_mer));
                                }
                                //System.out.println(Long.toBinaryString(k_mer));
                                addK_mers(kmersHash,k_mer,pins.get(i-1),L);
                            }

                        }
                    }
                }
                if (pins.size() < NUM) {
                    break;
                }
            } else {
                break;
            }
            pos = (Long) res.get("pos");
        }
        return kmersHash;
    }
    public Map<Long,K_mers> readFa(String filePath){
        final String ENCODING = "UTF-8";
        final int NUM = 4000;
        Map<Long,K_mers> kmersHash = new HashMap();  //k-mer s库
        File file = new File(filePath);
        long pos = 0L;
        String M ="A";
        int nu =0;
        while (true) {
            Map<String, Object> res = FileRead.BufferedRandomAccessFileReadLine(file, ENCODING, pos, NUM);
            // 如果返回结果为空结束循环
            if (org.apache.commons.collections4.MapUtils.isEmpty(res)) {
                break;
            }
            List<String> pins = (List<String>) res.get("pins");
            if (org.apache.commons.collections4.CollectionUtils.isNotEmpty(pins)) {
                //generate kmers
                for(int i=1; i<pins.size();i+=2){
                    //System.out.println(pins.get(i).length());
                    if(pins.get(i).contains("N")){
                        continue;
                    }
                    int count = (pins.get(i).length()-pins.get(i).replace(M, "").length())/M.length();
                    if(count/pins.get(i).length() > 0.8){
                        continue;
                    }
                    System.out.printf("%d---->length:%d\n",nu++,pins.get(i).length());
                    String read_seq = pins.get(i);  //get the right read
                    long k_mer = 0;
                    if(k < read_seq.length()){
                        for(int L=0; L<read_seq.length()-k+1;L++){
                            if(L==0){
                                k_mer = kmersAdders(read_seq.substring(0,k));
                                //System.out.println(Long.toBinaryString(k_mer));
                            }else{
                                k_mer =((k_mer <<2) | pre_base(read_seq.charAt(L+k-1)))&(addSum);
                                //System.out.println(Long.toBinaryString(k_mer));
                            }
                            //System.out.println(Long.toBinaryString(k_mer));
                            if(k_mer == 0L){
                                continue;
                            }
                            addK_mers(kmersHash,k_mer,pins.get(i-1),L);
                        }
                    }
                }
                if (pins.size() < NUM) {
                    break;
                }
            } else {
                break;
            }
            pos = (Long) res.get("pos");
        }
        return kmersHash;
    }
    public void addK_mers(Map<Long,K_mers>kmerhash,long base_seq,String id,int pos){
        boolean kmerin_read = false;

        if(kmerhash.containsKey(base_seq)){  //如果kmer表已经包含该kmer，则频率增加
            long fre = (long)kmerhash.get(base_seq).getFrequence();
            kmerhash.get(base_seq).setFrequence(fre+1);
            max_seq_fre(base_seq,fre+1);

            for(ReadID_Pos id_pos:kmerhash.get(base_seq).getReadID_Pos()){
                //如果当前kmer的readID_Pos记录中不包含该条read，则numOfRead+1
                if(id_pos.readID.equals(id)){
                    kmerin_read = true;
                    id_pos.readPos.add(pos);
                    break;
                }
            }
            if(!kmerin_read){
                kmerhash.get(base_seq).setNumOfRead(kmerhash.get(base_seq).getNumOfRead()+1);
                ReadID_Pos readID_pos = new ReadID_Pos();
                readID_pos.readID=id;
                readID_pos.readPos.add(pos);
                readID_pos.deleFlag = false;
                List<ReadID_Pos> list = kmerhash.get(base_seq).getReadID_Pos();
                //list.addAll(kmerhash.get(base_seq).getReadID_Pos());
                list.add(readID_pos);
                kmerhash.get(base_seq).setReadID_Pos(list);
                }
        }else{
            ReadID_Pos readID_pos = new ReadID_Pos();
            readID_pos.readID = id;
            readID_pos.readPos.add(pos);
            readID_pos.deleFlag = false;
            List<ReadID_Pos> list = new ArrayList<ReadID_Pos>();
            list.add(readID_pos);
            K_mers kmer = new K_mers(base_seq,list,1L,1L);
            kmerhash.put(base_seq,kmer);  //将所有kmers存入hash表
            max_seq_fre(base_seq,1);
        }
    }
    public K_mers createInitilKmer(Map<Long,K_mers> kmerhash){
        return kmerhash.get(seq_fre.seq);
    }
    public void merge_unipath(){
        Vertex start_vex = new Vertex();
        Vertex end_vex = new Vertex();
        int[] is_visit = new int[debru.vertexList.size()];
        for(int i=0; i<debru.vertexList.size();i++){
            Vertex vex = new Vertex();
            vex = debru.vertexList.get(i);
            //System.out.println(Integer.toString(i)+"-->in:"+Integer.toString(inDegree(vex))+"--out:"+Integer.toString(outDegree(vex)));
            if((is_visit[i]==0) && (debru.inDegree(vex)==1) && (debru.outDegree(vex)==1)){
                //is_visit[i] = 1;
                start_vex = vex;
                end_vex = vex;
                while (vex.firstOut != null){
                    vex = debru.vertexList.get((int)vex.firstOut.headvex);
                    if((is_visit[i]==0) && (debru.inDegree(vex)==1) && (debru.outDegree(vex)==1|debru.outDegree(vex)==0)){
                        is_visit[debru.vertexList.indexOf(vex)] =1;
                        if(debru.vertexList.get((int)vex.firstIn.tailvex) == start_vex){
                            end_vex = vex;
                            start_vex = debru.mergeVex_in_unipath(start_vex,end_vex);
                            end_vex = start_vex;
                            vex = start_vex;
                            is_visit[debru.vertexList.indexOf(end_vex)] =0;
                        }else{
                            start_vex = vex;
                            end_vex = vex;
                        }

                    }
                }
            }
        }
    }

    public long sumOf_2(int k){
        long sum =0L;
        for(int i=0;i<k;i++){
            sum = sum <<1 |1;
        }
        return sum;
    }
    public long pre_base(char c){
        long re =0;
        switch (c){
            case 'A':re =0;break;
            case 'C':re =1;break;
            case 'G':re =2;break;
            case 'T':re =3;break;
        }
        return re;
    }
    public long kmersAdders(String str){
        long adders = 0L;
        for(int i=0;i<str.length();i++){
            adders = adders <<2|pre_base(str.charAt(i));
        }

        return adders;
    }
    public void max_seq_fre(long seq,long fre){
        if(seq_fre.fre < fre){
            seq_fre.seq = seq;
            seq_fre.fre = fre;
        }
    }
    //.......///////
}