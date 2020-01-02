/**
 * 
 */
package ics.whu.edu.cn.madrix.clustering.wapper;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Administrator
 *
 */
public class ClusteringBenchTranslator {
    private double[][] data = null;
    private int[] labels;
    private int dim = -1;
    private int dis = -1;

    public ClusteringBenchTranslator(String path) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(path));
        String line = "";
        List<String> cache = new ArrayList<>();
        while ((line = br.readLine()) != null) {
            cache.add(line);
            /*String[] fields = line.split("\t");
            if (data == null) {
                dim = fields.length - 1;
                data = new double[dim][];
            }*/
        }
        br.close();
        data = new double[cache.size()][];
        labels = new int[cache.size()];
        int idx = 0;
        for (String str : cache) {
            String[] fields = str.split("\t");
            data[idx] = new double[fields.length + dis];
            for (int i = 0; i < fields.length + dis; i++) {
                data[idx][i] = Double.parseDouble(fields[i]);
            }
            labels[idx] = Integer.parseInt(fields[fields.length - 1]);
            idx++;
        }
        dim = data[0].length;
    }

    public double[][] getData() {
        return data;
    }

    public int[] getLabels() {
        return labels;
    }

    public int getDim() {
        return dim;
    }
}
