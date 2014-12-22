package org.bioinfo.babelomics.utils;

import org.bioinfo.commons.log.Logger;

public class TicTac {
	
	//LOCAL VARIABLES///////////////////////////////
	public static final short NO_VERBOSE = 0;
	public static final short VERBOSE = 1;
	public static final short VERBOSE_BEGIN = 2; 
	public static final short VERBOSE_END = 3;
    String name = "default";
    String initMessage = "begin of"; 
    String endMessage = "end of";
    long t1=0,t2=0,total=0;
    private short verbose = 3;
    boolean running = false;
    String els = "";
    private Logger out;
        
    //CONSTRUCTORS//////////////////////////////////
    public TicTac(){
    	//
	}
    public TicTac(String name){
		this.name=name;		
	}
    public TicTac(short verbose){
		this.verbose = verbose;		
	}
    public TicTac(Logger out){
		this.out=out;		
	}
	public TicTac(String name, short verbose){
		this.name=name;
		this.verbose = verbose;
	}	
	public TicTac(String name, Logger out){
		this.name=name;
		this.out=out;	
	}
	public TicTac(short verbose, Logger out){		
		this.verbose = verbose;
		this.out=out;		
	}
	public TicTac(String name, short verbose, Logger out){
		this.name=name;
		this.verbose = verbose;
		this.out=out;		
	}	
	
	//FUNCTIONS/////////////////////////////////////
	//Tic
	public void tic(){
		t1=System.currentTimeMillis();
		running = true;
		if ((verbose==VERBOSE)|(verbose==VERBOSE_BEGIN)){			
			print(">> " + initMessage + " '" + name + "'");			
		}
	}
	//Tic
	public void tic(String Name){
		name = Name;
		t1=System.currentTimeMillis();
		running = true;
		if ((verbose==VERBOSE)|(verbose==VERBOSE_BEGIN)){
			print(">> " + initMessage + " '" + name + "'");			
		}
	}
	//Tic
	public void tic(String Name, short v){
		name = Name;
		verbose = v;
		t1=System.currentTimeMillis();
		running = true;
		if ((verbose==VERBOSE)|(verbose==VERBOSE_BEGIN)){
			print(">> " + initMessage + " '" + name + "'");
		}
	}
	
	//Tac
	public void tac(){
		t2=System.currentTimeMillis();
		total = t2-t1;
		running = false;
		if ((verbose==VERBOSE)|(verbose==VERBOSE_END)){
			print("<< "  + endMessage + " '" + name + "' in " + getElapsedString());
		}		
	}
	
	//GETTERS///////////////////////////////////////
	//Get elapsed time
	public long getElapsed(){
		return total;
	}
	//Get elapsed time in string format
	public String getElapsedString(){
		long h,m,s,mm;
		double aux;
		String sh = "", sm = "", ss = "", smm = "";
		mm = total;
		if (mm>999){			
			aux = mm/1000;
			s = (long)Math.floor(aux);
			mm = mm - s*1000;
			smm = mm + "ms";
			if (s>59){
				aux = s/60;
				m = (long)Math.floor(aux);
				s = s - m*60;
				ss = s + "s ";
				if (m>59){
					aux = m/60;
					h = (long)Math.floor(aux);
					m = m - h*60;
					sm = m + "m ";
					sh = h + "h ";
				} else sm = m + "m ";
			}else ss = s + "s ";			
		} else smm = total + "ms";
		els = sh + sm + ss + smm;
		return els;
	}
	//Get status
	public boolean getStatus(){
		return running;
	}
	//Print output
	public void print(String message){
		if (out==null) System.out.println(message);
		else out.println(message);
	}
	//SETTERS///////////////////////////////////////
	//Set verbose
	public void setVerbose(short sv){
		verbose = sv;
	}
	//Set name
	public void setName(String Name){
		name = Name;
	}
	//Set total time
	public void setElapsed(long et){
		total = et;
	}
	public Logger getLogger() {
		return out;
	}
	public void setLogger(Logger out) {
		this.out = out;
	}
}
