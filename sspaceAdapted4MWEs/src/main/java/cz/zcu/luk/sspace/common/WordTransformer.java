package cz.zcu.luk.sspace.common;

import cz.zcu.luk.sspace.text.IteratorFactoryStopwords;

/**
 * Created with IntelliJ IDEA.
 * User: Lukr
 * Date: 16.2.13
 * Time: 15:07
 * To change this template use File | Settings | File Templates.
 */
public class WordTransformer {
    public static String getRealWord(String possibleWord) {
        String wordReal;
        if (possibleWord.endsWith(IteratorFactoryStopwords.STOPWORD_FLAG)) {
            wordReal = possibleWord.substring(0, possibleWord.length()-IteratorFactoryStopwords.STOPWORD_FLAG.length()); // store the stopword!
        }
        else {
            wordReal = possibleWord;
        }
        return wordReal;
    }

    public static String getWord(String possibleWord) {
        String word;
        if (possibleWord.endsWith(IteratorFactoryStopwords.STOPWORD_FLAG)) {
            word = IteratorFactoryStopwords.EMPTY_TOKEN;
        }
        else {
            word = possibleWord;
        }
        return word;
    }
}
