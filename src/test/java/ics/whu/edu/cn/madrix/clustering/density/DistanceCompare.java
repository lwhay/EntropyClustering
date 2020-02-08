package ics.whu.edu.cn.madrix.clustering.density;

import ics.whu.edu.cn.madrix.basis.MadrixUtils;
import ics.whu.edu.cn.madrix.common.exceptions.MadrixException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class DistanceCompare {
    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException, MadrixException {
        BufferedReader br = new BufferedReader(new FileReader("./src/test/java/ics/whu/edu/cn/madrix/clustering/resources/faceimage.txt"));
        double[] src = null, dest;
        int srclabel = 0, destlabel;
        String line;
        int number = 0;
        while ((line = br.readLine()) != null) {
            String[] fields = line.split("\t");
            if (number++ == 0) {
                src = new double[fields.length - 1];
                for (int i = 0; i < fields.length - 1; i++) {
                    src[i] = Double.parseDouble(fields[i]);
                }
                srclabel = Integer.parseInt(fields[fields.length - 1]);
            }

            dest = new double[fields.length - 1];
            for (int i = 0; i < fields.length - 1; i++) {
                dest[i] = Double.parseDouble(fields[i]);
            }
            destlabel = Integer.parseInt(fields[fields.length - 1]);
            System.out.println((number - 1) + ":" + MadrixUtils.vectorPNorm(src, dest, 2) + ":" + MadrixUtils.vectorPDistance(src, dest, 2) + ":" + MadrixUtils.vectorDTWDistance(src, dest) + "：" + srclabel + "-" + destlabel);
        }
        br.close();
    }
}
