package cn.contigs;

import org.graphstream.graph.Edge;
import org.graphstream.graph.Graph;
import org.graphstream.graph.implementations.MultiGraph;
import org.graphstream.ui.swingViewer.View;
import org.graphstream.ui.swingViewer.Viewer;
import scala.Int;

import java.util.*;

//import org.graphstream.algorithm.Toolkit

public class DeBruijnGraph {  //使用十字链表存储结构
    //以Kmer的hash表做表头，以Vertex做节点
    List<Vertex> vertexList = new ArrayList<>();
    Map<EdgeNode,Long> re_edgeNode = new HashMap<>();
    //int num_vertex = 0;
    public boolean addVex(K_mers kmer){
        if(getPosition(kmer,vertexList)<0){
            Vertex verx = new Vertex(kmer,null,null);
            vertexList.add(verx);
        }
        return true;
    }

    public boolean addEdge(K_mers preKmer,K_mers folKmer){
        EdgeNode edge = new EdgeNode();
        //EdgeNode edge_in = new EdgeNode();
        long vi = getPosition(preKmer,vertexList);
        long vj = getPosition(folKmer,vertexList);
        EdgeNode ed = new EdgeNode();
        ed = vertexList.get((int)vi).firstOut;
        while(ed !=null){
            if(ed.tailvex ==vi && ed.headvex==vj){
                if(!ed.idpos_map.isEmpty()){
                    for(String id : ed.idpos_map.keySet()){
                        if(!ed.idpos_map.get(id).isEmpty()){
                            ed.idpos_map.get(id).remove(0);
                            if(ed.idpos_map.get(id).isEmpty()){
                                ed.idpos_map.remove(id);
                            }
                            if(re_edgeNode.containsKey(ed)){
                                re_edgeNode.put(ed,re_edgeNode.get(ed)+1);
                            }else{
                                re_edgeNode.put(ed,1L);
                            }
                            return true;
                        }
                    }
                }else{
                    return false;
                }
            }
            ed = ed.taillink;
        }
        edge.tailvex = vi;
        edge.headvex = vj;
        edge.taillink = vertexList.get((int)vi).firstOut;
        edge.headlink = vertexList.get((int)vj).firstIn;
        edge.idpos_map = adjID_pos(vertexList.get((int)vi).kmer,vertexList.get((int)vj).kmer);
        vertexList.get((int)vi).firstOut = edge;
        vertexList.get((int)vj).firstIn = edge;

        //edge_in.tailvex = vi;
        //edge_in.headvex = vj;


        return true;
    }
    public Map<String,List<Integer>>adjID_pos(K_mers kmer1,K_mers kmer2){
        //too much new object   jfjafljalklkjajg

        Map<String,List<Integer>> idpos_map = new HashMap<>();
        List<String> idlist1 = new ArrayList<>();
        List<String> idlist2 = new ArrayList<>();

        for(ReadID_Pos idpos:kmer1.getReadID_Pos())    idlist1.add(idpos.readID);
        for(ReadID_Pos idpos:kmer2.getReadID_Pos())    idlist2.add(idpos.readID);
        idlist1.retainAll(idlist2);
        List<Integer> pos1 = new ArrayList<>();
        List<Integer> pos2 = new ArrayList<>();
        for(String id: idlist1){
            int index1 = index_ReadID_Pos(id,kmer1.getReadID_Pos());
            int index2 = index_ReadID_Pos(id,kmer2.getReadID_Pos());
            pos1 = kmer1.getReadID_Pos().get(index1).readPos;
            pos2 = kmer2.getReadID_Pos().get(index2).readPos;
            List<Integer> pos_list = new ArrayList<>();
            pos_list = dis_k_Ele(pos1,pos2,1);
            if(!pos_list.isEmpty()){
                idpos_map.put(id,pos_list);  // set the thtreahold is 1

            }
            pos1.clear();
            pos2.clear();
        }
        return idpos_map;
    }
    public List<Integer> dis_k_Ele(List<Integer> list1,List<Integer> list2,int h){
        int i=0 ,j=0;
        List<Integer> pos = new ArrayList<>();
        while(i<list1.size() && j<list2.size()){
            //System.out.println(i);
            if(list2.get((int)j) > list1.get((int)i) && list2.get((int)j) - list1.get((int)i) <= h){
                pos.add(list1.get((int)i));
                i++;
                j++;
            }
            else if(list1.get((int)i) < list2.get((int)j)) i++;
            else j++;
        }
        return pos;
    }
    public int index_ReadID_Pos(String id,List<ReadID_Pos> id_posList){
        for(int i=0; i<id_posList.size();i++){
            if(id_posList.get(i).readID.equals(id)){
                return i;
            }
        }
        return -1;
    }
    public Vertex mergeVex_in_unipath(Vertex vi,Vertex vj){
        Vertex v = new Vertex();
        K_mers kmer = new K_mers();
        long vi_fre = vi.kmer.getFrequence();
        long vj_fre = vj.kmer.getFrequence();
        kmer.setFrequence(vi_fre>vj_fre?vj_fre:vi_fre);
        long vi_numRead = vi.kmer.getNumOfRead();
        long vj_numRead = vj.kmer.getNumOfRead();
        kmer.setNumOfRead(vi_numRead>vj_numRead?vj_numRead:vi_numRead);
        long vi_seq = vi.kmer.getbase_seq();
        long vj_seq = vj.kmer.getbase_seq();
        kmer.setbase_seq((vi_seq << 2)|(vj_seq & 3));

        List<String> id_list = new ArrayList<>();
        Map<String,List<Integer>> map = new HashMap<>();
        for(ReadID_Pos r_p :vi.kmer.getReadID_Pos()){
            map.put(r_p.readID,r_p.readPos);
        }
        for(ReadID_Pos r_p :vj.kmer.getReadID_Pos()){
            if(map.containsKey(r_p.readID)){
                id_list.add(r_p.readID);
            }
        }
        List read_posList = new ArrayList();
        for(String id :id_list){
            ReadID_Pos read_pos = new ReadID_Pos();
            read_pos.readID = id;
            read_pos.deleFlag = false;
            read_pos.readPos = map.get(id);
            read_posList.add(read_pos);
        }
        kmer.setReadID_Pos(read_posList);
        v.kmer = kmer;
        v.firstIn = vi.firstIn;
        v.firstOut = vj.firstOut;
        if(v.firstOut!=null){
            v.firstOut.tailvex = vi.firstOut.tailvex;
        }
        re_edgeNode.put(v.firstOut,re_edgeNode.get(vi.firstOut));
        re_edgeNode.remove(vi.firstOut);
        re_edgeNode.remove(vj.firstOut);
        vj.firstIn = null;
        vj.firstOut = null;
        //vertexList.get((int)vj.firstOut.headvex).firstIn.headlink = v.firstOut;
        vertexList.set((int)vertexList.indexOf(vi),v);

        return v;
    }
    public long getPosition(K_mers kmer,List<Vertex> verx){
        for(Vertex v :verx){
            if(v.kmer.equals(kmer)){
                return verx.indexOf(v);
            }
        }
        return -1;
    }
    public int outDegree(Vertex v){
        int out_count =0;
        int index = vertexList.indexOf(v);
        if(vertexList.get(index).firstOut !=null){
            out_count++;
            EdgeNode mEdgeNode = new EdgeNode();
            mEdgeNode = vertexList.get(index).firstOut;
            while (mEdgeNode.taillink != null) {
                mEdgeNode = mEdgeNode.taillink;
                out_count++;
            }
        }
        return out_count;
    }
    public int inDegree(Vertex v){
        int in_count =0;
        int index = vertexList.indexOf(v);
        if(vertexList.get(index).firstIn !=null){
            EdgeNode mEdgeNode = new EdgeNode();
            mEdgeNode = vertexList.get(index).firstIn;
            in_count++;
            while (mEdgeNode.headlink!= null) {
                mEdgeNode = mEdgeNode.headlink;
                in_count++;
            }
        }
        return in_count;
    }
    public void showGraph(String id){
        Graph graph = new MultiGraph(id);
        for(Vertex v :vertexList){
            if(v.firstOut==null && v.firstIn ==null) continue;
            graph.addNode(Long.toString(vertexList.indexOf(v)));
        }
        for(int i=0;i<vertexList.size();i++){
            if(vertexList.get(i).firstOut !=null){

                EdgeNode mEdgeNode = new EdgeNode();
                mEdgeNode = vertexList.get(i).firstOut;
                String id_vi = Long.toString(mEdgeNode.tailvex);
                String id_vj = Long.toString(mEdgeNode.headvex);
//                String vi = Long.toString(vertexList.indexOf(vertexList.get((int)edgeNode.tailvex)));
//                String vj = Long.toString(vertexList.indexOf(vertexList.get((int)edgeNode.headvex)));
                graph.addEdge(id_vi+id_vj, id_vi, id_vj,true);
                while (mEdgeNode.taillink != null) {
                    mEdgeNode = mEdgeNode.taillink;
                    id_vi = Long.toString(mEdgeNode.tailvex);
                    id_vj = Long.toString(mEdgeNode.headvex);
//                    vi = Long.toString(vertexList.indexOf(vertexList.get((int)edgeNode.tailvex)));
//                    vj = Long.toString(vertexList.indexOf(vertexList.get((int)edgeNode.headvex)));
                    graph.addEdge(id_vi+id_vj, id_vi, id_vj,true);
                }
            }
        }
        graph.addAttribute("ui.stylesheet","url(./style/stylesheet)");
        graph.addAttribute("ui.quality");
        graph.addAttribute("ui.antialias");
        System.setProperty("gs.ui.renderer", "org.graphstream.ui.j2dviewer.J2DGraphRenderer");
        Viewer viewer = graph.display();
        View view = viewer.getDefaultView();
        //view.resizeFrame(8000, 6000);
        //view.setViewCenter(440000, 2503000, 0);
        //view.setViewPercent(0.85);

    }
    public void printf_reEdge(){
        for(Map.Entry<EdgeNode,Long>entry : re_edgeNode.entrySet()){
            System.out.println(Long.toString(entry.getKey().tailvex) +"-------->"+
                    Long.toString(entry.getKey().headvex)+":  "+Long.toString(entry.getValue()));
        }
    }
    public void printf(){
        for(int i=0;i<vertexList.size();i++){
            System.out.print(Integer.toString(i)+"-->");
            if(vertexList.get(i).firstOut !=null){
                EdgeNode mEdgeNode = new EdgeNode();
                mEdgeNode = vertexList.get(i).firstOut;
                System.out.print(Long.toString(mEdgeNode.headvex)+"--");
                while (mEdgeNode.taillink != null) {
                    mEdgeNode = mEdgeNode.taillink;
                    System.out.print(Long.toString(mEdgeNode.headvex)+"--");
                }
            }
            System.out.println();
        }

    }

}
class Vertex{
    K_mers kmer;
    EdgeNode firstIn;
    EdgeNode firstOut;

    public Vertex() { }

    public Vertex(K_mers kmer, EdgeNode firstIn, EdgeNode firstOut) {
        this.kmer = kmer;
        this.firstIn = firstIn;
        this.firstOut = firstOut;
    }
}
class EdgeNode{
    long tailvex;
    long headvex;
    EdgeNode headlink;
    EdgeNode taillink;
    Map<String,List<Integer>> idpos_map = new HashMap<>();
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        EdgeNode edgeNode = (EdgeNode) o;
        return tailvex == edgeNode.tailvex &&
                headvex == edgeNode.headvex &&
                Objects.equals(headlink, edgeNode.headlink) &&
                Objects.equals(taillink, edgeNode.taillink);
    }

    @Override
    public int hashCode() {

        return Objects.hash(tailvex, headvex, headlink, taillink);
    }
}