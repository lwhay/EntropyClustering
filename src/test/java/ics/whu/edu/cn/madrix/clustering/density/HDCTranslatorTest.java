/**
 * 
 */
package ics.whu.edu.cn.madrix.clustering.density;

import java.io.IOException;

import ics.whu.edu.cn.madrix.clustering.wapper.HDCTranslator;

/**
 * @author Administrator
 *
 */
public class HDCTranslatorTest {

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        HDCTranslator hdct = new HDCTranslator(args[0], args[1]);
        hdct.dummy();
    }

}
