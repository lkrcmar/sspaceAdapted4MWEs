package cz.zcu.luk.sspace.common;

import java.io.*;

/**
 * Created with IntelliJ IDEA.
 * User: Lukr
 * Date: 13.4.13
 * Time: 11:53
 * To change this template use File | Settings | File Templates.
 */
public class Serializer {

    public static void serializeData(Object data, String outputFN) {

        FileOutputStream fos;
        try {
            fos = new FileOutputStream(outputFN);
            ObjectOutput oos = new ObjectOutputStream(fos);
            oos.writeObject(data);
            oos.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println("Serialized: " + outputFN);
    }

    public static Object deserialiazeData(String outputFN) {
        // if file does not exist, return null..
        if (!(new File(outputFN).exists())) return null;

        FileInputStream fis;
        Object data = null;
        try {
            fis = new FileInputStream(outputFN);
            ObjectInputStream ois = new ObjectInputStream(fis);
            data = ois.readObject();
            ois.close();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }

        System.out.println("Deserialized: " + outputFN);
        return data;
    }
}
