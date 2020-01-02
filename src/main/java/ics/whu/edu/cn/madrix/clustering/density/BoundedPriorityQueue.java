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

public class BoundedPriorityQueue<E> extends PriorityQueue<E> {
    private static final long serialVersionUID = 1L;
    private int limit;

    public BoundedPriorityQueue(int maxCapacity, boolean asc) {
        super(maxCapacity, asc ? null : Collections.reverseOrder());
        this.limit = maxCapacity;
    }

    public E pop() {
        E r = null;
        int cur = 0;
        while (cur++ < limit && (r = super.poll()) != null) {

        }
        return r;
    }

    public int getLimit() {
        return limit;
    }
}