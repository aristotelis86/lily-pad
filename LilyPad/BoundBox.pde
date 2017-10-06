/**********************************************************************
      BoundBox class: Creates bounding boxes for the springs.
      They are used in collision detection. Both AABB and OBB
      boxes can be created.

Example code:
      (To be filled)

**********************************************************************/
class BoundBox {
  Spring spA;
  LineSegment [] lines = new LineSegment[4];
  PVector [] vertices = new PVector[5];
  Window view;
  PVector Min, Max;
  OrthoNormal orth;
  PVector Normal;
  
  
  BoundBox( Spring sp1_ ) {
    spA = sp1_;
    this.BasicInfo();
    this.createAABB();
    view = sp1_.myWindow;
  }
  
  void BasicInfo() {
    orth = new OrthoNormal( spA.p1.position, spA.p2.position );
    Normal = new PVector(orth.nx, orth.ny);
    Normal.setMag(spA.thick/2.);
    
    float xc = 0.5*(spA.p1.position.x + spA.p2.position.x);
    float yc = 0.5*(spA.p1.position.y + spA.p2.position.y);
    LineSegment line1 = new LineSegment( PVector.add(spA.p1.position, Normal), PVector.add(spA.p2.position, Normal));
    LineSegment line2 = new LineSegment( PVector.sub(spA.p1.position, Normal), PVector.sub(spA.p2.position, Normal));
    
    lines[0] = line1;
    lines[1] = new LineSegment( new PVector(line1.Start.x, line1.Start.y), new PVector(line2.Start.x, line2.Start.y));
    lines[2] = line2;
    lines[3] = new LineSegment( new PVector(line2.End.x, line2.End.y), new PVector(line1.End.x, line1.End.y));
    vertices[0] = new PVector(xc, yc);
    vertices[1] = new PVector(line1.Start.x, line1.Start.y);
    vertices[2] = new PVector(line1.End.x, line1.End.y);
    vertices[3] = new PVector(line2.Start.x, line2.Start.y);
    vertices[4] = new PVector(line2.End.x, line2.End.y);
    
  }
  
  void createAABB() {  
    float [] Xs = {lines[0].Start.x, lines[0].End.x, lines[2].Start.x, lines[2].End.x};
    float [] Ys = {lines[0].Start.y, lines[0].End.y, lines[2].Start.y, lines[2].End.y};
    
    Min = new PVector();
    Min.x = min(Xs);
    Min.y = min(Ys);
    Max = new PVector();
    Max.x = max(Xs);
    Max.y = max(Ys);
  }
  
  void displayOBB() {
    stroke(255, 0, 0);
    for (LineSegment ls : lines) {
      line(view.px(ls.Start.x), view.py(ls.Start.y), view.px(ls.End.x), view.py(ls.End.y));
    }
  }
  
  void displayAABB() {
    noStroke();
    fill(0, 255, 0);
    rect(view.px(Min.x), view.py(Min.y), view.px(Max.x-Min.x), view.py(Max.y-Min.y));
  }
}

/**********************************************************************
      LineSegment class: Creates a line segment given its
      starting and ending points (PVectors). Auxiliary 
      class to the collision detection methods.

Example code:
      (To be filled)

**********************************************************************/
class LineSegment {
  PVector Start;
  PVector End;
  OrthoNormal orth;
  
  LineSegment( PVector s, PVector e ) {
    Start = s;
    End = e;
    orth = new OrthoNormal( Start, End );
  }
}