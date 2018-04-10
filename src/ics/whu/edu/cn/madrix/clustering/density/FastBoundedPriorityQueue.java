/*
 * Copyright 2007-2015 by The Regents of the Wuhan University of China.
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * you may obtain a copy of the License from
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package ics.whu.edu.cn.madrix.clustering.density;

import java.util.Collections;
import java.util.PriorityQueue;

public class FastBoundedPriorityQueue<E> extends PriorityQueue<E> {
    private static final long serialVersionUID = 1L;
    private int limit;
    private E peek;
    private boolean asc = true;

    public FastBoundedPriorityQueue(int maxCapacity, boolean asc) {
        super(maxCapacity, asc ? Collections.reverseOrder() : null);
        this.limit = maxCapacity;
        this.asc = asc;
    }

    @SuppressWarnings("unchecked")
    @Override
    public boolean add(E e) {
        if (super.size() < limit) {
            boolean ret = super.add(e);
            peek = super.peek();
            return ret;
        } else {
            if (asc) {
                if (((Comparable<E>) peek).compareTo(e) > 0) {
                    super.remove();
                    boolean ret = super.add(e);
                    peek = super.peek();
                    return ret;
                }
            } else {
                if (((Comparable<E>) peek).compareTo(e) < 0) {
                    super.remove();
                    boolean ret = super.add(e);
                    peek = super.peek();
                    return ret;
                }
            }
        }
        return false;
    }

    public int getLimit() {
        return limit;
    }
}