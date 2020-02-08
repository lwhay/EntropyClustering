/**
 * 
 */
package ics.whu.edu.cn.madrix.clustering.density;

import ics.whu.edu.cn.madrix.common.exceptions.MadrixException;

/**
 * @author Administrator
 *
 */
public interface IClustering {
    public int[] export();

    public void action() throws MadrixException;
}
