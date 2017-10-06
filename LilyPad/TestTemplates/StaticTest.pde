int nx = 150; // x-dir resolution
int ny = 150; // y-dir resolution

float Length = 20;
float thick = 1;
float MassNum = 0.1;
int resol = 1;
float stiffness = 100;
PVector lpos = new PVector(nx/2.,4);
PVector align = new PVector(0,1);
FlexibleSheet sheet;
Window view; // convert pixels to non-dim frame
WriteInfo writer;

PVector gravity = new PVector(0,5);
float t = 0;
float dt;

void settings(){
    size(600, 600);
}

void setup() {
  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
  sheet = new FlexibleSheet( Length, thick, MassNum, resol, stiffness, lpos, align, view );
  sheet.cpoints[0].makeFixed();
  
  sheet.Calculate_Stretched_Positions( gravity );
  
  dt = sheet.dtmax;
  
  writer = new WriteInfo(sheet);
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
