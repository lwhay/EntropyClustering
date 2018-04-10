/*
 * Copyright 2007-2016 by The Regents of the Wuhan University of China.
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
package ics.whu.edu.cn.madrix.common.exceptions;

/**
 * @author Michael
 *
 */
public class MadrixException extends Exception {

    private static final long serialVersionUID = 1L;

    public MadrixException() {
    }

    /**
     * @param message
     */
    public MadrixException(String message) {
        super(message);
    }

    /**
     * @param cause
     */
    public MadrixException(Throwable cause) {
        super(cause);
    }

    /**
     * @param message
     * @param cause
     */
    public MadrixException(String message, Throwable cause) {
        super(message, cause);
    }

    /**
     * @param message
     * @param cause
     * @param enableSuppression
     * @param writableStackTrace
     */
    public MadrixException(String message, Throwable cause, boolean enableSuppression, boolean writableStackTrace) {
        super(message, cause, enableSuppression, writableStackTrace);
    }

}