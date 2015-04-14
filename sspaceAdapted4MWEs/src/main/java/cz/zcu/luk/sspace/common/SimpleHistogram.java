package cz.zcu.luk.sspace.common;

/**
 * Created with IntelliJ IDEA.
 * User: Lukr
 * Date: 20.1.13
 * Time: 18:15
 * To change this template use File | Settings | File Templates.
 */
public class SimpleHistogram {
    // in zero position the amount of values > maxValue is stored
    private int[] freqCounts;
    private int maxValue;

    public SimpleHistogram(int maxValue) {
        freqCounts = new int[maxValue+1];
        this.maxValue = maxValue;
    }

    public void addOccurrenceOfValue(int value) {
        if (value > maxValue) {
            freqCounts[0]++;
        }
        else {
            freqCounts[value]++;
        }
    }

    public void printHistogram() {
        for (int i = 1; i < freqCounts.length; i++) {
            System.out.println(i + "\t" + freqCounts[i]);
        }
        System.out.println(">" + maxValue + "\t" + freqCounts[0]);
    }
}
