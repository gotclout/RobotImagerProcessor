package image.processor;

import imagehunter.connection.HunterSender;
import imagehunter.server.ByteStream;

import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.Socket;

import javax.imageio.ImageIO;

public class RobotConnection {
	
	private boolean traveling = false;
	private Socket commandSocket;
	private Socket dataSocket;
	private Socket cameraSocket;
	private BufferedReader dataReader;
	private PrintWriter commandSender;
	
	
	public RobotConnection(String host){
		try{
			commandSocket = new Socket(host, 2699);
			dataSocket = new Socket(host, 2700);
			cameraSocket = new Socket(host, 2750);
			dataReader = new BufferedReader(new InputStreamReader(dataSocket.getInputStream()));
			commandSender = new PrintWriter(commandSocket.getOutputStream(), true);
			DataThread dt = new DataThread();
			dt.start();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}

	public void goTo(double x, double y, double angle){
		commandSender.println("goto " + x + " " + y + " " + angle);
		commandSender.flush();
	}
	
	public void pan(int degrees){
		commandSender.println("p " + degrees);
		commandSender.flush();
	}
	
	public void tilt(int degrees){
		commandSender.println("i " + degrees);
		commandSender.flush();
	}
	
	public BufferedImage getImg(){
		BufferedImage ret = null;
		try{
			File f = new File("IMG.png");
			if(f.exists()){
				f.delete();
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
		try{
			ByteStream.toStream(cameraSocket.getOutputStream(), "TAKEPIC");
			File f = new File("IMG.png");
			//sleep here and give the agent time to save and write the image
			//to the stream how long?
			Thread.sleep(10000);
			ByteStream.toFile(cameraSocket.getInputStream(), f);
			ret = ImageIO.read(f);
		}
		catch(Exception e){
			e.printStackTrace();
		}
		return ret;
	}
	
	public void stop(){
		
	}
	
	public boolean isTraveling(){
		return traveling;
	}
	
	private class DataThread extends Thread{
		public void run(){
			while(true){
				try{
				    while (!dataReader.ready()) {}
					String line = dataReader.readLine();
					if(line.contains("true")){
						traveling = true;
					}
					else{
						traveling = false;
					}
				}
				catch(Exception e){
					e.printStackTrace();
				}
			}
		}
	}
	
	public static void main(String args[]){
		RobotConnection con = new RobotConnection("192.168.0.100");
		HunterSender sender = new HunterSender("192.168.7.126");
		BufferedImage img = con.getImg();
		sender.sendImageToServer(0, 0, 0.0f, "test1", img);
		System.out.println("Width: " + img.getWidth());
		img = con.getImg();
		sender.sendImageToServer(0, 0, 0.0f, "test2", img);
		System.out.println("Width: " + img.getWidth());
		img = con.getImg();
		sender.sendImageToServer(0, 0, 0.0f, "test3", img);
		System.out.println("Width: " + img.getWidth());
		try{
			ImageIO.write(img, "png", new File("book_test.png"));
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
}
