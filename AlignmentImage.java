/*
 * IMPORTANT: arguments from the command line should be filenames - INPUT and OUTPUT respectively
 */

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import javax.imageio.ImageIO;


/**
 *
 * @author Tristan
 */
//T - green
//A - red
//G - purple
//C - yellow
//indel - white (blank)
public class AlignmentImage {

    public static int width = 0;
    public static int height = 1000;
    public static final int nuc_width = 1;
    public static final int nuc_height = 1;
    
    static ArrayList<String> list = new ArrayList<>();

    public static void fillRect(BufferedImage im, int startX, int startY, int w, int h, Color c) {
        for (int i = startX; i < startX + w; i++) {
            for (int j = startY; j < startY + h; j++) {
                im.setRGB(i, j, c.getRGB());
            }
        }
    }

    public static void fillImage(BufferedImage im, ArrayList<String> strings) {
//        int a_count = 0, c_count = 0, g_count = 0, t_count = 0;

        Graphics2D graphics = im.createGraphics();
        graphics.setColor(Color.WHITE);
        graphics.fillRect(0, 0, im.getWidth(), im.getHeight());

        for (int i = 0; i < strings.size(); i++) {
            for (int j = 0; j < strings.get(i).length(); j++) {
                if (strings.get(i).substring(j, j + 1).equals("A")) {
                    fillRect(im, 50 + nuc_width * j, 50 + nuc_height * i, nuc_width, nuc_height, Color.RED);
                    //a_count++;
                } else if (strings.get(i).substring(j, j + 1).equals("G")) {
                    fillRect(im, 50 + nuc_width * j, 50 + nuc_height * i, nuc_width, nuc_height, new Color(190, 0, 95));
                    //g_count++;
                } else if (strings.get(i).substring(j, j + 1).equals("T")) {
                    fillRect(im, 50 + nuc_width * j, 50 + nuc_height * i, nuc_width, nuc_height, Color.GREEN);
                    //t_count++;
                } else if (strings.get(i).substring(j, j + 1).equals("C")) {
                    fillRect(im, 50 + nuc_width * j, 50 + nuc_height * i, nuc_width, nuc_height, new Color(255, 255, 0));
                    //c_count++;
                }
//                fillRect(im,50+nuc_width*j, 50+nuc_height*i, nuc_width, 1, Color.BLACK);
            }
        }
//        System.out.println(a_count+" "+c_count+" "+g_count+" "+t_count);
//        System.out.println((double) a_count/(a_count+c_count+g_count+t_count));
    }

    public static void main(String[] args) throws IOException {

        File input = new File(args[0]);
        BufferedReader reader = new BufferedReader(new FileReader(input));

        String currentLine;
        String currentSeq = "";

        currentLine = reader.readLine();
        while (currentLine != null) {
            if(currentLine.length()==0){
                
            }
            else if (currentLine.substring(0, 1).equals(">")) {
                if (currentSeq.length()>0) {
                    list.add(currentSeq);
                    if(currentSeq.length()*nuc_width+50>=width)
                        width = currentSeq.length()*nuc_width+100;
                }
                currentSeq = "";
            } else {
                currentSeq += currentLine;
            }
            currentLine = reader.readLine();
        }
        
        list.add(currentSeq);
        if(currentSeq.length()*nuc_width+50>=width){
            width = currentSeq.length()*nuc_width+100;
        }
        height = 100 + list.size()*nuc_height;
        
        BufferedImage img = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

        fillImage(img, list);

//Makes BitMap file
//        try {
//            ImageIO.write(img, "BMP", new File("ouput.bmp"));
//        } catch (Exception IOException) {
//            System.out.println("IOException caught.");
//        }
        
        try {
            ImageIO.write(img, "PNG", new File(args[1]));
        } catch (IOException e) {
            System.out.println("Exception caught.");
        }

    }

}
