package cz.zcu.luk.sspace.config;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Lukr
 * Date: 10.10.12
 * Time: 8:21
 * To change this template use File | Settings | File Templates.
 */
public class Config {
    private static Config instance = null;

    private final String configFN = System.getProperty("configFN");
    // param X value
    public final Map<String, String> configuration;

    protected Config() {
        configuration = new HashMap<String, String>();
        String line = null;
        try {
            BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(configFN), "UTF-8"));
            while ((line = br.readLine()) != null) {
                if (!line.startsWith("#") && !line.isEmpty()) {
                    configuration.put(line.trim().split("=")[0], line.trim().split("=")[1]);
                }
            }
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static Config getInstance() {
        if (instance == null) {
            instance = new Config();
        }
        return instance;
    }
}
