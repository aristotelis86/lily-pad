int nx = 150; // x-dir resolution
int ny = 150; // y-dir resolution

float Length = ny/6.;
float thick = 1;
float MassNum = 0.5;
int resol = 1;
float stiffness = 300;
PVector lpos = new PVector(nx/2.,ny/20.);
PVector align = new PVector(0,1);
FlexibleSheet [] sheet;
Window view; // convert pixels to non-dim frame
WriteInfo writer;
CollisionHandler collide;

PVector gravity = new PVector(0,10);
float t = 0;
float dt;

void settings(){
    size(600, 600);
}

void setup() {
  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
  sheet = new FlexibleSheet[2];
  
  sheet[0] = new FlexibleSheet( Length, thick, MassNum, resol, stiffness, lpos, align, view );
  sheet[0].cpoints[0].makeFixed();
  sheet[0].Calculate_Stretched_Positions( gravity );
  
  sheet[1] = new FlexibleSheet( nx/2., thick, 4*MassNum, resol, 4*stiffness, new PVector(nx/4.,0.55*ny), new PVector(1,0), view );
  sheet[1].cpoints[0].makeFixed();
  sheet[1].cpoints[sheet[1].numOfpoints-1].makeFixed();
  
  dt = sheet[0].dtmax;
  if (sheet[1].dtmax<dt) dt = sheet[1].dtmax;
  
  writer = new WriteInfo(sheet);
  for (int i=0; i<2; i++) {
    writer.addGenInfo( i, "Gravity is "+gravity.mag()+" in y-direction");
    writer.addGenInfo( i, "Stretched length is "+ sheet[i].CurrentLength() +" in y-direction");
    writer.addGenInfo( i, "Stretching ratio is "+ (sheet[i].CurrentLength()-sheet[i].Length)/sheet[i].Length);
    writer.addGenInfo( i, "Flexible sheets catch test.");
  }
  
  collide = new CollisionHandler( sheet );
} // end of setup

void draw() {
  background(185); 
  fill(0); // color of text for timer
  textSize(32); // text size of timer
  text(t, 10, 30); // position of timer
  
  for (FlexibleSheet fs : sheet) {
    fs.updateAlt( dt, gravity );
    fs.updateAlt2( dt, gravity );
  }
  
  collide.HandleCollisions(); 
  
  for (FlexibleSheet fs : sheet) fs.display(color(0, 255, 0));
  
  writer.dampAllInfo( t, true );
  
  if (t>5) sheet[0].cpoints[0].makeFree();
  t += dt;
}


void keyPressed() {
  writer.terminateFiles();
}


