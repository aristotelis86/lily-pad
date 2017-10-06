int nx = 150; // x-dir resolution
int ny = 150; // y-dir resolution

float Length = ny/4.;
float thick = 1;
float MassNum = 0.5;
int resol = 1;
float stiffness = 500;
PVector lpos = new PVector(nx/3.,ny/2.);
PVector align = new PVector(0,1);
FlexibleSheet sheet;
Window view; // convert pixels to non-dim frame
WriteInfo writer;

PVector gravity = new PVector(0,10);
float t = 0;
float dt;

float maxVel = 20;

void settings(){
    size(600, 600);
}

void setup() {
  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
  sheet = new FlexibleSheet( Length, thick, MassNum, resol, stiffness, lpos, align, view );
  sheet.cpoints[0].makeFixed();
  sheet.Calculate_Stretched_Positions( gravity );
  
  // Apply the impulse
  int N = sheet.numOfpoints;
  float [] vx = new float[N];
  float [] vy = new float[N];
  
  for (int i = 1; i < N; i++) {
    vx[i] = ((i-1)/(N-2)) * maxVel;
    vy[i] = sheet.cpoints[i].velocity.y;
    sheet.cpoints[i].UpdateVelocity( vx[i], vy[i] );
  }
  
  dt = sheet.dtmax;
  
  writer = new WriteInfo(sheet);
  writer.addGenInfo( 0, "Gravity is "+gravity.mag()+" in y-direction");
  writer.addGenInfo( 0, "Stretched length is "+ sheet.CurrentLength() +" in y-direction");
  writer.addGenInfo( 0, "Stretching ratio is "+ (sheet.CurrentLength()-sheet.Length)/sheet.Length);
  writer.addGenInfo( 0, "Velocity of trailing edge is "+ maxVel +" in x-direction.");
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


