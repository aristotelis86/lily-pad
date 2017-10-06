int nx = 150; // x-dir resolution
int ny = 150; // y-dir resolution

float Length = ny/4.;
float thick = 1;
float MassNum = 0.5;
int resol = 1;
float stiffness = 100;
PVector lpos = new PVector(nx/3.,ny/2.);
PVector align = new PVector(0,1);
FlexibleSheet sheet;
Window view; // convert pixels to non-dim frame
WriteInfo writer;

PVector gravity = new PVector(0,0);
float t = 0;
float dt;

float sinN = 1; // mode

void settings(){
    size(600, 600);
}

void setup() {
  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
  sheet = new FlexibleSheet( Length, thick, MassNum, resol, stiffness, lpos, align, view );
  sheet.cpoints[0].makeFixed();
  sheet.cpoints[sheet.numOfpoints-1].makeFixed();
  
  // Create the distortion
  float L = sheet.Length;
  float sinAmp = L/15; // amplitude of init disturbance
  int N = sheet.numOfpoints;
  float [] x = new float[N];
  float [] y = new float[N];
  x[0] = lpos.x;
  y[0] = lpos.y;
  for (int i = 1; i < N; i++) {
    y[i] = (L/(N-1)) + y[i-1];
    x[i] = sinAmp * sin(sinN*PI*(y[i]-lpos.y)/L) + lpos.x;
  }
  sheet.UpdateState( x, y );
  
  dt = sheet.dtmax;
  
  writer = new WriteInfo(sheet);
  writer.addGenInfo( 0, "Gravity is "+gravity.mag()+" in y-direction");
  writer.addGenInfo( 0, "Stretched length is "+ sheet.CurrentLength() +" in y-direction");
  writer.addGenInfo( 0, "Stretching ratio is "+ (sheet.CurrentLength()-sheet.Length)/sheet.Length);
  writer.addGenInfo( 0, "Normal mode "+ sinN +" is tested.");
} // end of setup

void draw() {
  background(185); 
  fill(0); // color of text for timer
  textSize(32); // text size of timer
  text(t, 10, 30); // position of timer
  
  sheet.updateAlt( dt, gravity );
  sheet.updateAlt2( dt, gravity );
  sheet.display(color(0, 255, 0));
  
  writer.dampAllInfo( t, true );
  
  t += dt;
}


void keyPressed() {
  writer.terminateFiles();
}


