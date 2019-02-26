package cn.contigs;

import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class test {
    public static void main(String args[]) {
        List<Integer> list  = new ArrayList();
        Scanner sc = new Scanner(System.in);
        int input = 1;
        while (input != 0){
            input = sc.nextInt();
            if(input != 0){
                list.add(input);
            }
        }
        for(int i=0;i<list.size();i++){
            System.out.println(numofWater(list.get(i),3));
        }

    }
    public static int numofWater(int n,int m){
        int num = 0;
        while(n/m >0){
            num += n/m;
            n = n/m + n%m;
        }
        if(n%m ==2){
            num ++;
        }
        return num;
    }
}
