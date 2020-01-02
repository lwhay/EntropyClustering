/**
 * 
 */
package ics.whu.edu.cn.madrix.clustering.wapper;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Administrator
 *
 */
public class HDCTranslator {

    public HDCTranslator(String path, String opath) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(path));
        String line = "";
        List<String> cache = new ArrayList<>();
        while ((line = br.readLine()) != null) {
            cache.add(line);
        }
        br.close();
        Map<String, Integer> labels = new HashMap<>();
        BufferedWriter bw = new BufferedWriter(new FileWriter(opath));
        int count = 0;
        for (String str : cache) {
            String os = "";
            String[] fields = str.split(",");
            for (int i = 0; i < fields.length - 1; i++) {
                os += fields[i];
                os += "\t";
            }
            String label = fields[fields.length - 1];
            if (labels.containsKey(label)) {
                os += labels.get(label);
            } else {
                os += count;
                labels.put(label, count++);
            }
            os += "\n";
            bw.write(os);
        }
        bw.close();
    }

    public void dummy() {

    }
}
