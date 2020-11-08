picture axis(real x0,real x1,pen e=black+linewidth(0.25mm),
             pen t=black+linewidth(0.1mm),real divs=1,int step=1,real scaled=1.0,
             real unitsize=1,real scaletext=0.667,real scaleticks=0.5,int digits=4) {
  picture p;	
  real fx = scaleticks*1.5/((x1-x0)*unitsize^2); 
  real unit = 1.0/divs; 
  pair r = (x0-unit/2.0,0); 
  pair s = (x1+unit/2.0,0);
  draw(p,r--s,e);
  real unit = 1/divs; 
  int steps = 0;
  for(real i=x0;i<=x1;i+=unit) {
  	++steps;
    int signe = 1;
    if(i > 0) signe = -1;
    draw(p,(i,0)--(i,-2*fx),e);
    if(steps % step == 0)  
      label(p,scale(scaletext)*("$"+string(scaled*i,digits)+"$"),(i,-3*fx),S,e);
    if ((i != x0) || (i < 0)) 
      for(int k=1;k<=3;++k)
        draw(p,(i+signe*0.25*k*unit,0)--(i+signe*0.25*k*unit,-fx),t);
  }
  return p;
}

picture axes(real x0,real x1,real y0,real y1,pen e=black+linewidth(0.25mm),
              pen t=black+linewidth(0.1mm),real divs=1,int step=1,
              real unitsize=1,real scaletext=0.667,real scaleticks=0.5,int digits=4) {
  picture p;	
  real fx = scaleticks*min(1.5/((x1-x0)*unitsize^2),1.5/((y1-y0)*unitsize^2)); 
  path r = (x0,y0)--(x1,y0)--(x1,y1)--(x0,y1)--cycle; 
  draw(p,r,e);
  real unit = 1/divs; 
  int steps = 0;
  for(real i=x0;i<=x1;i+=unit) {
  	++steps;
    int signe = 1;
    if(i > 0) signe = -1;
    draw(p,(i,y0)--(i,y0-2*fx),e);
    draw(p,(i,y1)--(i,y1+2*fx),e);
    if(steps % step == 0) {
      real ii = i;
      if(fabs(i) < 0.001) ii = 0.0; 
      label(p,scale(scaletext)*("$"+string(ii,digits)+"$"),(i,y0-3*fx),S,e);
      label(p,scale(scaletext)*("$"+string(ii,digits)+"$"),(i,y1+3*fx),N,e);
    }
    if ((i != x0) || (i < 0)) 
      for(int k=1;k<=3;++k) {
        draw(p,(i+signe*0.25*k*unit,y0)--(i+signe*0.25*k*unit,y0-fx),t);
        draw(p,(i+signe*0.25*k*unit,y1)--(i+signe*0.25*k*unit,y1+fx),t);
    }
  }
  steps = 0;
  for(real i=y0;i<=y1;i+=unit) {
  	++steps;
  	int signe = 1;
    if(i > 0) 
      signe = -1;
    draw(p,(x0,i)--(x0-2*fx,i),e);
    draw(p,(x1,i)--(x1+2*fx,i),e);
    if(steps % step == 0) { 
      real ii = i;
      if(fabs(i) < 0.001) ii = 0.0; 
      label(p,scale(scaletext)*("$"+string(ii,digits)+"$"),(x0-3*fx,i),W,e);
      label(p,scale(scaletext)*("$"+string(ii,digits)+"$"),(x1+3*fx,i),E,e);
    }
    for(int k=1;k<=3;++k) {
      draw(p,(x0,i+signe*0.25*k*unit)--(x0-fx,i+signe*0.25*k*unit),t);
      draw(p,(x1,i+signe*0.25*k*unit)--(x1+fx,i+signe*0.25*k*unit),t);
    }
  }
  return p;
}

picture canonica(real x0,real x1,real y0,real y1,pen e=black+linewidth(0.25mm),
              pen t=black+linewidth(0.1mm),real divs=1,int step=1,
			 real unitsize=1,real scaletext=0.667,real scaleticks=0.5,int digits=4,bool base=false) {
  picture p;
  real fx = scaleticks*min(1.5/((x1-x0)*unitsize^2),1.5/((y1-y0)*unitsize^2));
  pair r = (x0,0); 
  pair s = (x1,0);
  draw(p,r--s,e);
  r = (0,y0); 
  s = (0,y1);
  draw(p,r--s,e);

  real unit = 1/divs; 
  int steps = 0;
  for(real i=x0;i<=x1;i+=unit) {
  	++steps;
    int signe = 1;
    if(i > 0) signe = -1;
    draw(p,(i,0)--(i,-2*fx),e);
    if(steps % step == 0) {
      real ii = i;
      if(fabs(i) < 0.001) ii = 0.0; 
      if (ii != 0)
        label(p,scale(scaletext)*("$"+string(ii,digits)+"$"),(i,-3*fx),S,e);
      if ((i != x0) || (i < 0)) 
        for(int k=1;k<=3;++k) {
          draw(p,(i+signe*0.25*k*unit,0)--(i+signe*0.25*k*unit,-fx),t);
        }
    }
  }
  steps = 0;
  for(real i=y0;i<=y1;i+=unit) {
  	++steps;
  	int signe = 1;
    if(i > 0) 
      signe = -1;
    draw(p,(0,i)--(-2*fx,i),e);
    if(steps % step == 0) { 
      real ii = i;
      if(fabs(i) < 0.001) ii = 0.0;
      if (ii != 0)
        label(p,scale(scaletext)*("$"+string(ii,digits)+"$"),(-3*fx,i),W,e);
      for(int k=1;k<=3;++k) {
        draw(p,(0,i+signe*0.25*k*unit)--(-fx,i+signe*0.25*k*unit),t);
      }
    }
  }
  if(base) {
	r = (0,0); 
	s = (1,0);
    draw(p,r--s,e+linewidth(0.4mm),Arrow(size=1.3mm));
    s = (0,1);
	draw(p,r--s,e+linewidth(0.4mm),Arrow(size=1.3mm));
  }
  return p;
}

picture grid(real x0,real x1,real y0,real y1,pen e=gray(0.7)+linewidth(0.1mm),
             pen t=gray(0.55)+linewidth(0.067mm),int divs=1) {
  picture p;
  real unit = 1/divs;
	for(real i=y0;i<=y1;i+=unit) {
    if(i < y1) {
      for(int k=1;k<=3;++k) {
        if (i+0.25*k*unit <= y1)
         draw(p,(x0,i+0.25*k*unit)--(x1,i+0.25*k*unit),t);
      }
    }
  }
  for(real i=x0;i<=x1;i+=unit) {
    if(i < x1) {
      for(int k=1;k<=3;++k) {
        draw(p,(i+0.25*k*unit,y0)--(i+0.25*k*unit,y1),t);
      }
    }
  }
	for(real i=y0;i<=y1;i+=unit) {
    draw(p,(x0,i)--(x1,i),e);
  }
	for(real i=x0;i<=x1;i+=unit) {
	  draw(p,(i,y0)--(i,y1),e);
  }
  return p;
}

void square(picture p,pair center,real size,pen color,real degrees=0) {
  pair A, B, C, D;
  A = center - (size/2,size/2);
  B = center + (size/2,-size/2);
  C = center + (size/2,size/2);
  D = center + (-size/2,size/2);
  filldraw(p,rotate(degrees,center)*(A--B--C--D--cycle),color+opacity(0.2),color+2bp);
}
