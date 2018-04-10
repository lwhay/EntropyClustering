package ics.whu.edu.cn.madrix.clustering.density;

import java.util.Comparator;
import java.util.Map;
import java.util.Map.Entry;

class MapValueComparator implements Comparator<Map.Entry<Integer,Double>> {

    @Override
    public int compare(Entry<Integer,Double> me1, Entry<Integer,Double> me2) {

        return me1.getValue().compareTo(me2.getValue());
    }
}

