package org.bactosim.haldane;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javax.vecmath.Point2i;

import org.apache.commons.math3.geometry.euclidean.threed.Line;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import java.awt.geom.Point2D;
import java.awt.geom.Point2D.Double;
import com.sun.security.auth.UnixNumericUserPrincipal;

import repast.simphony.context.Context;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.continuous.NdPoint;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.util.ContextUtils;
import repast.simphony.util.collections.IndexedIterable;
import repast.simphony.visualization.visualization2D.Random2DLayout;

@SuppressWarnings({"rawtypes", "unchecked"})
public class EnginePhy {
	private double x;
	private double y;
	private int xmax;
	private int ymax;

	private Grid grid;
	private ContinuousSpace space;
	private Class clazz = null;

	public EnginePhy(double xx, double yy, int xxmax, int yymax) {
		x = xx;
		y = yy;
		xmax = xxmax;
		ymax = yymax;
		clazz = (new VEColi()).getClass(); 
	}

	//@ScheduledMethod(start = 1, interval = 2, shuffle=true)
	public void step1() {
		Context context = (Context) ContextUtils.getContext(this);
		IndexedIterable<VEColi> agents = context.getObjects(clazz);

		space = (ContinuousSpace) context.getProjection("continuous-space");
		grid = (Grid) context.getProjection("grid-space");

		GridPoint p;
		NdPoint point;

		//for (Object o: grid.getObjectsAt(p.getX(),p.getY())){
		for (VEColi agent : agents) {
			if(agent.getElongation() > 0.6) {
				p = grid.getLocation(agent);
				point = space.getLocation(agent);distance(p.getX(), p.getY());
			}
		}
	}

	@ScheduledMethod(start = 1, interval = 10, shuffle=true)
	public void step() {
		RingVisitor();
	}

	public void RingVisitor() {
		int counter= 0;
		double rings = Math.min(Math.min(xmax - x, x), Math.min(ymax - y, y));
		ColonyRing cr= new ColonyRing((int) x, (int)y);

		// Ring loop
		for(int ring = 0; ring < rings; ring++) {
			// Ring components loop

			counter= ColonyEdge.agentCount(this, cr.Items(ring));
			if(ring > 1 && counter == 0) break;

			ColonyExpansion.Relaxation(this, cr.Items(ring));
		}


	}

	public int count(int x, int y) {
		int v= 0;
		for (@SuppressWarnings("unused") Object o: grid.getObjectsAt(x, y)) {
			v++;
		}
		return v;
	}

	public double distance(int x, int y) {
		double v= 0;
		if(count(x,y) > 1) {
			for (Object a1: grid.getObjectsAt(x, y)) {
				for (Object a2: grid.getObjectsAt(x, y)) {
					if(!a1.equals(a2)) {
						v= Agent.distance((VEColi) a1, (VEColi) a2);
						
					}
				}
			}

		}
		return v;
	}

}


/**
 * This class calculates the distance between two 
 * agents
 * TODO: Improve and make it a public class     
 */
@SuppressWarnings({"rawtypes"})
class Agent {

	public static double heading(Vector3D v) {
		return Angle.convheading(Math.toDegrees(Math.atan2((v.getY()),(v.getX()))));
	}

	public static double heading0(VEColi a1) {
		Vector3D v1 = vector0(a1);
		Vector3D v2 = vector1(a1); 
		return heading(v2.subtract(v1).normalize());
	}

	public static Vector3D vector(VEColi a1) {
		Vector3D point0 = vector0(a1);
		Vector3D point1 = vector1(a1); 
		return point1.subtract(point0);
	}

	public static double heading11(VEColi a1) {
		Line l = line(a1);

		double alpha = Math.atan2( l.getDirection().getY() , l.getDirection().getX());
		double alphadeg = Math.toDegrees(alpha);
		double h = Angle.convheading(alphadeg);
		return h;
	}

	public static Line line(VEColi a1) {
		Vector3D v1 = vector0(a1);
		Vector3D v2 = vector1(a1);
		return (new Line(v1, v2, 1.0e-10));
	}


	public static Vector3D vector0(VEColi a1) {
		NdPoint p;
		Context context = (Context) ContextUtils.getContext(a1);
		ContinuousSpace space = (ContinuousSpace) context.getProjection("continuous-space");
		p= space.getLocation(a1);
		Vector3D v1 = new Vector3D(p.getX(), p.getY(), 0);
		return v1;
	}

	public static Vector3D vector1(VEColi a1) {
		Vector3D p = vector0(a1);
		double alpha = Math.toRadians(Angle.convheading(a1.getHeading()));
		Vector3D v2 = new Vector3D(p.getX() + a1.getLength() *  Math.cos(alpha), p.getY() + a1.getLength() * Math.sin(alpha), 0);
		return v2;
	}

	public static boolean collided(VEColi a1, VEColi a2, double overlap) {
		double closest = a1.getWidth() * (1 - overlap);
		return (distance(a1, a2) < closest);
	}



	public static double distance(VEColi a1, VEColi a2) {
		Vector3D uS = vector0(a1);
		Vector3D uE = vector1(a1);
		Vector3D vS = vector0(a2);
		Vector3D vE = vector1(a2);
		Vector3D w1 = vector0(a1);
		Vector3D w2 = vector0(a2);

		Vector3D u = uE.subtract(uS);
		Vector3D v = vE.subtract(vS);
		Vector3D w = w1.subtract(w2);

		double a = u.dotProduct(u);
		double b = u.dotProduct(v);
		double c = v.dotProduct(v);
		double d = u.dotProduct(w);
		double e = v.dotProduct(w);
		double D = a * c - b * b;
		double sc, sN, sD = D;
		double tc, tN, tD = D;

		if (D < 0.01) {
			sN = 0;
			sD = 1;
			tN = e;
			tD = c;
		} else {
			sN = (b * e - c * d);
			tN = (a * e - b * d);
			if (sN < 0) {
				sN = 0;
				tN = e;
				tD = c;
			} else if (sN > sD) {
				sN = sD;
				tN = e + b;
				tD = c;
			}
		}

		if (tN < 0) {
			tN = 0;
			if (-d < 0) {
				sN = 0;
			} else if (-d > a) {
				sN = sD;
			} else {
				sN = -d;
				sD = a;
			}
		} else if (tN > tD) {
			tN = tD;
			if ((-d + b) < 0) {
				sN = 0;
			} else if ((-d + b) > a) {
				sN = sD;
			}
			else {
				sN = (-d + b);
				sD = a;
			}
		}

		if (Math.abs(sN) < 0.01) {
			sc = 0;
		} else {
			sc = sN / sD;
		}

		if (Math.abs(tN) < 0.01) {
			tc = 0;
		} else {
			tc = tN / tD;
		}

		Vector3D dP = w.add(u.scalarMultiply(sc).subtract((v.scalarMultiply(tc))));
		return Math.sqrt(dP.dotProduct(dP));
	}


}


class Angle {

	public static double convheading(double alpha) {
		return (450.0 - alpha) % 360.0;
	}
}





/**
 * This class provides Shoving on ring  
 * TODO: Improve and make it a public class     
 */
class ColonyExpansion {


	public static Point2D.Double interseccion(Point2D.Double a ,Point2D.Double b, Point2D.Double c, Point2D.Double d) {
		Point2D.Double P1 = a;
		Point2D.Double P2 = b;
		Point2D.Double P3 = c;
		Point2D.Double P4 = d;
		double x1 = P1.getX();
		double y1 = P1.getY();
		double x2 = P2.getX();
		double y2 = P2.getY();
		double x3 = P3.getX();
		double y3 = P3.getY();
		double x4 = P4.getX();
		double y4 = P4.getY();

		double m1 = (y2 - y1) / ( x2 - x1);
		double n1 = y1 - m1 * x1;
		double m2 = (y4 - y3) / ( x4 - x3);
		double n2 = y3 - m2 * x3;
		double A = -m1 ;
		double B = 1;
		double C = n1;
		double D = -m2;
		double E = 1;
		double F = n2;
		double x,y;


		if((x2-x1) == 0) {
			x = x1;
			if((y4-y3) == 0) {
				y = y4;
			} else {
				y = m2 * x + n2;
			}

		} else if ((y2 - y1) == 0 ) {
			y = y1;
			if((x4-x3) == 0 ) {
				x = x4;
			} else {
				x = (y-n2) / m2;
			}
		} else {
			y = ((F*A) - (D*C)) / ((E*A) - (D*B));
			x = (C - (B*y)) / A;
		}



		Point2D.Double sol = new Point2D.Double(x, y);
		return sol;
	}

	//Para calcular el area utilizaremos la formula de Heron.
	public static double calculateArea (Point2D.Double A, Point2D.Double B, Point2D.Double C) {

		Point2D.Double P1 = A; 
		Point2D.Double P2 = B; 
		Point2D.Double P3 = C;

		double a = Math.sqrt( Math.pow((P2.getX() - P1.getX()),2) + Math.pow((P2.getY() - P1.getY()),2) ); //lado que une A y B
		double b = Math.sqrt( Math.pow((P3.getX() - P2.getX()),2) + Math.pow((P3.getY() - P2.getY()),2) ); //lado que une B y C
		double c = Math.sqrt( Math.pow((P3.getX() - P1.getX()),2) + Math.pow((P3.getY() - P1.getY()),2) ); //lado que une A y C

		double area = Math.sqrt( (a+b+c) * (a+b-c) * (b+c-a) * (c+a-b) ) / 4;

		return area;
	}




	//Metodo para obtener las coordenadas del baricentro de un triangulo.
	public static Point2D.Double getBar (Point2D.Double p1, Point2D.Double p2 , Point2D.Double p3) {
		Point2D.Double A = p1;
		Point2D.Double B = p2;
		Point2D.Double C = p3;

		return new Point2D.Double((A.getX() + B.getX() + C.getX())/3 , (A.getY() + B.getY() + C.getY())/3);

	}



	//static double h = 0;
	@SuppressWarnings("rawtypes")
	public static void Relaxation(Object o, Iterable<Point2i> ring) {
		Context context = (Context) ContextUtils.getContext(o);
		Grid grid = (Grid) context.getProjection("grid-space");
		ContinuousSpace space = (ContinuousSpace) context.getProjection("continuous-space");

		// Loop on ring boxes
		for(Point2i p: ring) {

			for (Object a1: grid.getObjectsAt(p.getX(), p.getY())) {
				for (Object a2: grid.getObjectsAt(p.getX(), p.getY())) {
					if(!(a1 instanceof VEColi) || !(a2 instanceof VEColi) || a1.equals(a2)) continue;
					if( Agent.collided( (VEColi) a1, (VEColi) a2, 0.1) ) {

						//Teniendo en cuenta que la bacteria a1 es el centro de coordenadas tiene los puntos A1, B1, C1, D1
						double xA1 = - (( (VEColi) a1).getLength())/2;  
						double yA1 = (( (VEColi) a1).getWidth())/2;
						double xB1 = - (( (VEColi) a1).getLength())/2;
						double yB1 = - (( (VEColi) a1).getWidth())/2;
						double xC1 = (( (VEColi) a1).getLength())/2;
						double yC1 = - (( (VEColi) a1).getWidth())/2;
						double xD1 = (( (VEColi) a1).getLength())/2;
						double yD1 = (( (VEColi) a1).getWidth())/2;


						//Angulo entre bacteria a1 y bacteria a2
						double alpha1 = Math.toRadians(Angle.convheading(((VEColi) a1).getHeading())); //angulo de a1 en radianes
						double alpha2 = Math.toRadians(Angle.convheading(((VEColi) a2).getHeading())); //angulo de a2 en radianes
						double alpha = alpha1 - alpha2;
						
						NdPoint a1Location = space.getLocation(a1); //Posicion de la celula a1.
						double xa1 = a1Location.getX();
						double ya1 = a1Location.getY();
						NdPoint a2Location = space.getLocation(a2); //Posicion de la celula a2.
						double xa2 = a2Location.getX();
						double ya2 = a2Location.getY();

						//double xG2 = Math.cos(alpha) * (xa2 - xa1) + Math.sin(alpha) * (ya2 - ya1);
						//Dado que el centro de G1 con el cambio de coordenadas esta en (0,0).
						double xG2 = Math.cos(alpha) * xa2 + Math.sin(alpha) * ya2;
						double yG2 = -Math.sin(alpha) * xa2 + Math.cos(alpha) * ya2;

						// Puntos de la Bacteria a2 A2,B2,C2,D2 teniendo en cuenta que el origen de coordenadas esta en el centro de la Bacteria1.
						double xA2 = xG2 - ( ((VEColi) a2).getLength()/2) * Math.cos( alpha ) - ( (( (VEColi) a2).getWidth())/2 * Math.sin( alpha ));
						double yA2 = yG2 - ( ((VEColi) a2).getLength()/2) * Math.sin( alpha ) + ( (( (VEColi) a2).getWidth())/2 * Math.cos( alpha ));

						double xB2 = xG2 - ( ((VEColi) a2).getLength()/2) * Math.cos( alpha ) + ( (( (VEColi) a2).getWidth())/2 * Math.sin( alpha ));
						double yB2 = yG2 - ( ((VEColi) a2).getLength()/2) * Math.sin( alpha ) - ( (( (VEColi) a2).getWidth())/2 * Math.cos( alpha ));

						double xC2 = xG2 + ( ((VEColi) a2).getLength()/2) * Math.cos( alpha ) + ( (( (VEColi) a2).getWidth())/2 * Math.sin( alpha ));
						double yC2 = yG2 + ( ((VEColi) a2).getLength()/2) * Math.sin( alpha ) - ( (( (VEColi) a2).getWidth())/2 * Math.cos( alpha ));

						double xD2 = xG2 + ( ((VEColi) a2).getLength()/2) * Math.cos( alpha ) - ( (( (VEColi) a2).getWidth())/2 * Math.sin( alpha ));
						double yD2 = yG2 + ( ((VEColi) a2).getLength()/2) * Math.sin( alpha ) + ( (( (VEColi) a2).getWidth())/2 * Math.cos( alpha ));

						//Voy a ver donde se encuentra la bacteria 2 con respecto a la bacteria 1 para poder hacer el cambio.
						boolean cambio = false;
						int caso = 0;
						if( xG2 > 0 && yG2 > 0 && alpha > 0 ){
							cambio = false;
							caso = 1;
						}else if( xG2 > 0 && yG2 > 0 && alpha < 0 ) {
							//Caso 2.  Simetria respecto al EJE X. 
							yA2 = -yA2;
							yB2 = -yB2;
							yC2 = -yC2;
							yD2 = -yD2;
							cambio = true;
							caso = 2; //Me encuentro en el caso2.
						} else if(xG2 < 0 && yG2 > 0 && alpha > 0) {
							cambio = false;
							caso = 3;
						}else if( xG2 < 0 && yG2 > 0 && alpha < 0 ) {
							//Caso 4. Simetria respecto al EJE Y.
							xA2 = -xA2;
							xB2 = -xB2;
							xC2 = -xC2;
							xD2 = -xD2;
							cambio = true;
							caso = 4; //Me encuentro en el caso4.
						} else if( xG2 < 0 && yG2 < 0 && alpha > 0 ) {
							//Caso 5. Simetria respecto al EJE X y al EJE Y.
							xA2 = -xA2;
							xB2 = -xB2;
							xC2 = -xC2;
							xD2 = -xD2;
							yA2 = -yA2;
							yB2 = -yB2;
							yC2 = -yC2;
							yD2 = -yD2;
							cambio = true;
							caso = 5;
						} else if( xG2 < 0 && yG2 < 0 && alpha < 0 ) {
							//Caso 6. Simetria respecto al EJE X.
							yA2 = -yA2;
							yB2 = -yB2;
							yC2 = -yC2;
							yD2 = -yD2;
							cambio = true;
							caso = 6; //Me encuentro en el caso5.

						} else if( xG2 > 0 && yG2 < 0 && alpha > 0 ){
							cambio = false;
							caso = 7;
						} else if( xG2 > 0 && yG2 < 0 && alpha < 0 ) {
							//Caso 8. Simetria respecto al EJE Y.
							xA2 = -xA2;
							xB2 = -xB2;
							xC2 = -xC2;
							xD2 = -xD2;
							cambio = true;
							caso = 8; //Me encuentro en el caso8.
						}



						Point2D.Double A1 = new Point2D.Double(xA1, yA1);
						Point2D.Double B1 = new Point2D.Double(xB1, yB1);
						Point2D.Double C1 = new Point2D.Double(xC1, yC1);
						Point2D.Double D1 = new Point2D.Double(xD1, yD1);
						Point2D.Double A2 = new Point2D.Double(xA2, yA2);
						Point2D.Double B2 = new Point2D.Double(xB2, yB2);
						Point2D.Double C2 = new Point2D.Double(xC2, yC2);
						Point2D.Double D2 = new Point2D.Double(xD2, yD2);


						Point2D.Double P = interseccion(A1,D1,A2,B2);
						Point2D.Double Q = interseccion(A1,D1,B2,C2);
						Point2D.Double R = interseccion(C1,D1,B2,C2);
						Point2D.Double S = interseccion(C1,D1,A2,B2);
						Point2D.Double T = interseccion(A1,B1,A2,B2);
						Point2D.Double U = interseccion(A1,B1,B2,C2);
						Point2D.Double V = interseccion(A1,D1,A2,D2);
						Point2D.Double W = interseccion(A1,B1,A2,D2);
						Point2D.Double X = interseccion(B1,C1,A2,B2);
						Point2D.Double Y = interseccion(B1,C1,B2,C2);
						Point2D.Double Z = interseccion(C1,D1,A2,D2);
						Point2D.Double N = interseccion(B1,C1,A2,D2);


						double solapamiento = 0.01 ;
						Point2D.Double baricentro = new Point2D.Double(0,0);

						if( yA2>yA1 && yB1<yB2 && yB2<yA1) {
							if( P.getX() < xA1) {
								if(xB2 < xA1) {
									//*** Caso 1 ***
									//Calculo el area de UA1Q
									solapamiento = solapamiento + calculateArea(U,A1,Q);
									//Calculo el baricentro
									baricentro = getBar(U,A1,Q);
								} else if(xB2 > xA1) {
									//*** Caso 2 ***
									//Calculo el area de TA1B2 + A1QB2
									double areaTA1B2 = calculateArea(T,A1,B2);
									double areaA1QB2 = calculateArea(A1,Q,B2);
									solapamiento = solapamiento + areaTA1B2 + areaA1QB2;
									//Calculo el baricentro
									Point2D.Double bar1 = getBar(T,A1,B2);
									Point2D.Double bar2 = getBar(A1,Q,B2);
									baricentro.setLocation((solapamiento * (bar1.getX() + bar2.getX()))/solapamiento, (solapamiento * (bar1.getY() + bar2.getY()))/solapamiento);
								}
							} else if(P.getX() > xA1) {
								if(xB2 > xD1) {
									//*** Caso 3 ***
									//Calculo el area de PD1S
									solapamiento = solapamiento + calculateArea(P,D1,S);
									//Calculo el baricentro
									baricentro = getBar(P,D1,S);
								} else if(xB2 < xD1) {
									if( Q.getX() < xD1 ) {
										//*** Caso 4 ***
										//Calculo el area de PB2Q
										solapamiento = solapamiento + calculateArea(P,B2,Q);
										//Calculo el baricentro
										baricentro = getBar(P,B2,Q);
									} else if( Q.getX() > xD1 ) {
										//*** Caso 5 ***
										//Calculo el area de PB2D1 y B2RD1
										double areaPB2D1 = solapamiento + calculateArea(P,B2,D1);
										double areaB2RD1 = solapamiento + calculateArea(B2,R,D1);
										solapamiento = areaPB2D1 + areaB2RD1;
										//Calculo el baricentro
										Point2D.Double bar1 = getBar(P,B2,D1);
										Point2D.Double bar2 = getBar(B2,R,D1);
										baricentro.setLocation((solapamiento * (bar1.getX() + bar2.getX()))/solapamiento, (solapamiento * (bar1.getY() + bar2.getY()))/solapamiento);
									}
								}
							}
						} else if( yB1<yA2 && yA2 < yA1 && yB1<yB2 && yB2<yA1 ){
							if(xC1 > xB2) {
								if(xA2<xA1) {
									//*** Caso 6 ***
									//Calcular el area de TWV + TVB2 + VB2Q
									double areaTWV = calculateArea(T,W,V);
									double areaTVB2 = calculateArea(T,V,B2);
									double areaVB2Q = calculateArea(V,B2,Q);
									solapamiento = solapamiento + areaTWV + areaTVB2 + areaVB2Q;
									//Calculo el baricentro
									Point2D.Double bar1 = getBar(T,W,V);
									Point2D.Double bar2 = getBar(T,V,B2);
									Point2D.Double bar3 = getBar(V,B2,Q);
									baricentro.setLocation((solapamiento * (bar1.getX()+bar2.getX()+bar3.getX()))/solapamiento, (solapamiento * (bar1.getY()+bar2.getY()+bar3.getY()))/solapamiento);
								} else if(xA2>xA1) {
									if(Q.getX() < xD1) {
										//*** Caso 7 ***
										//Calcular el area de A2VB2 + VB2Q
										double areaA2VB2 = calculateArea(A2,V,B2);
										double areaVB2Q = calculateArea(V,B2,Q);
										solapamiento = solapamiento + areaA2VB2 + areaVB2Q;
										//Calcular el baricentro
										Point2D.Double bar1 = getBar(A2,V,B2);
										Point2D.Double bar2 = getBar(V,B2,Q);
										baricentro.setLocation((solapamiento * (bar1.getX() + bar2.getX()))/solapamiento, (solapamiento * (bar1.getY() + bar2.getY()))/solapamiento);
									} else if(Q.getX() > xD1) {
										//*** Caso 8 ***
										//Calcular el area de VA2B2 + B2VD1 + D1RB2
										double areaVA2B2 = calculateArea(V,A2,B2);
										double areaB2VD1 = calculateArea(B2,V,D1);
										double areaD1RB2 = calculateArea(D1,R,B2);
										solapamiento = solapamiento + areaVA2B2 + areaB2VD1 + areaD1RB2;

										//Calcular el baricentro
										Point2D.Double bar1 = getBar(V,A2,B2);
										Point2D.Double bar2 = getBar(B2,V,D1);
										Point2D.Double bar3 = getBar(D1,R,B2);
										baricentro.setLocation((solapamiento * (bar1.getX()+bar2.getX()+bar3.getX()))/solapamiento, (solapamiento * (bar1.getY()+bar2.getY()+bar3.getY()))/solapamiento);
									}
								}
							} else if(xC1 < xB2) {
								if( V.getX() < xD1 ) {
									//*** Caso 9 ***
									//Calcular el area de A2VD1 + A2SD1
									double areaA2VD1 = calculateArea(A2,V,D1);
									double areaA2SD1 = calculateArea(A2,S,D1);
									solapamiento = solapamiento + areaA2VD1 + areaA2SD1;
									//Calcular el baricentro
									Point2D.Double bar1 = getBar(A2,V,D1);
									Point2D.Double bar2 = getBar(A2,S,D1);
									baricentro.setLocation((solapamiento * (bar1.getX() + bar2.getX()))/solapamiento, (solapamiento * (bar1.getY() + bar2.getY()))/solapamiento);

								} else if( V.getX() > xD1 ) {
									//*** Caso 10 ***
									//Calcular el area de ZA2S
									solapamiento = solapamiento + calculateArea(Z,A2,S);
									//Calcular el baricentro
									baricentro = getBar(Z,A2,S);
								}
							}
						} else if( yB1<yA2 && yA2<yA1 && yB2<yB1) {
							if(xA2 > xA1) {
								if(Y.getX()<xC1) {
									if(Q.getX() < xD1) {
										//*** Caso 11 ***
										//Calcular el area de XA2V + XVQ + XYQ
										double areaXA2V = calculateArea(X,A2,V);
										double areaXVQ = calculateArea(X,V,Q);
										double areaXYQ = calculateArea(X,Y,Q);
										solapamiento = solapamiento + areaXA2V + areaXVQ + areaXYQ;
										//Calcular el baricentro
										Point2D.Double bar1 = getBar(X,A2,V);
										Point2D.Double bar2 = getBar(X,V,Q);
										Point2D.Double bar3 = getBar(X,Y,Q);
										baricentro.setLocation((solapamiento * (bar1.getX()+bar2.getX()+bar3.getX()))/solapamiento, (solapamiento * (bar1.getY()+bar2.getY()+bar3.getY()))/solapamiento);

									} else if(Q.getX() > xD1) {
										//*** Caso 12 ***
										//Calcular el area de A2XV + XVY + VD1Y + YD1R
										double areaA2XV = calculateArea(A2,X,V);
										double areaXVY = calculateArea(X,V,Y);
										double areaVD1Y = calculateArea(V,D1,Y);
										double areaYD1R = calculateArea(Y,D1,R);
										solapamiento = solapamiento + areaA2XV + areaXVY + areaVD1Y + areaYD1R;
										//Calcular el baricentro
										Point2D.Double bar1 = getBar(A2,X,V);
										Point2D.Double bar2 = getBar(X,V,Y);
										Point2D.Double bar3 = getBar(V,D1,Y);
										Point2D.Double bar4 = getBar(Y,D1,R);
										baricentro.setLocation((solapamiento * (bar1.getX()+bar2.getX()+bar3.getX()+bar4.getX()))/solapamiento, (solapamiento * (bar1.getY()+bar2.getY()+bar3.getY()+bar4.getY()))/solapamiento);

									}
								} else if(Y.getX()>xC1) {
									if(X.getX() < xC1) {
										//*** Caso 13***
										//Calcular el area de VA2X + VXD1 + D1XC1
										double areaVA2X = calculateArea(V,A2,X);
										double areaVXD1 = calculateArea(V,X,D1);
										double areaD1XC1 = calculateArea(D1,X,C1);
										solapamiento = solapamiento + areaVA2X + areaVXD1 + areaD1XC1;
										//Calcular el baricentro
										Point2D.Double bar1 = getBar(V,A2,X);
										Point2D.Double bar2 = getBar(V,X,D1);
										Point2D.Double bar3 = getBar(D1,X,C1);
										baricentro.setLocation((solapamiento * (bar1.getX()+bar2.getX()+bar3.getX()))/solapamiento, (solapamiento * (bar1.getY()+bar2.getY()+bar3.getY()))/solapamiento);

									} else if(X.getX() > xC1) {
										//*** Caso 10 ***
										//Calcular el area de ZA2S
										solapamiento = solapamiento + calculateArea(Z,A2,S);
										//Calcular el baricentro
										baricentro = getBar(Z,A2,S);
									}
								}
							} else if(xA2 < xA1) {
								if(X.getX() < xB1) {
									//*** Caso 14 ***
									//Calcular el area de VWB1 + B1VY + VYQ
									double areaVWB1 = calculateArea(V,W,B1);
									double areaB1VY = calculateArea(B1,V,Y);
									double areaVYQ = calculateArea(V,Y,Q);
									solapamiento = solapamiento + areaVWB1 + areaB1VY + areaVYQ;
									//Calcular el baricentro
									Point2D.Double bar1 = getBar(V,W,B1);
									Point2D.Double bar2 = getBar(B1,V,Y);
									Point2D.Double bar3 = getBar(V,Y,Q);
									baricentro.setLocation((solapamiento * (bar1.getX()+bar2.getX()+bar3.getX()))/solapamiento, (solapamiento * (bar1.getY()+bar2.getY()+bar3.getY()))/solapamiento);


								} else if (X.getX() > xB1) {
									if(V.getX() < xA1) {
										//*** Caso 15 ***
										//Calcular el area de TA1Q + TQX + QXY
										double areaTA1Q = calculateArea(T,A1,Q);
										double areaTQX = calculateArea(T,Q,X);
										double areaQXY = calculateArea(Q,X,Y);
										solapamiento = solapamiento + areaTA1Q + areaTQX + areaQXY;
										//Calcular el baricentro
										Point2D.Double bar1 = getBar(T,A1,Q);
										Point2D.Double bar2 = getBar(T,Q,X);
										Point2D.Double bar3 = getBar(Q,X,Y);
										baricentro.setLocation((solapamiento * (bar1.getX()+bar2.getX()+bar3.getX()))/solapamiento, (solapamiento * (bar1.getY()+bar2.getY()+bar3.getY()))/solapamiento);

									} else if( V.getX() > xA1 ) {
										//*** Caso 16 ***
										//Calcular el area de XA2V + XVQ + XYQ
										double areaXA2V = calculateArea(X,A2,V);
										double areaXVQ = calculateArea(X,V,Q);
										double areaXYQ = calculateArea(X,Y,Q);
										solapamiento = solapamiento + areaXA2V + areaXVQ + areaXYQ;
										//Calcular el baricentro
										Point2D.Double bar1 = getBar(X,A2,V);
										Point2D.Double bar2 = getBar(X,V,Q);
										Point2D.Double bar3 = getBar(X,Y,Q);
										baricentro.setLocation((solapamiento * (bar1.getX()+bar2.getX()+bar3.getX()))/solapamiento, (solapamiento * (bar1.getY()+bar2.getY()+bar3.getY()))/solapamiento);

									}
								}
							}
						} else if( yA2<yB1 && yB2<yB1) {
							if(Y.getX() > xC1) {
								//*** Caso 17 ***
								//Calcular el area de NC1V + VC1D1
								double areaNC1V = calculateArea(N,C1,V);
								double areaVC1D1 = calculateArea(V,C1,D1);
								solapamiento = solapamiento + areaNC1V + areaVC1D1; 
								//Calcular el baricentro
								Point2D.Double bar1 = getBar(N,C1,V);
								Point2D.Double bar2 = getBar(V,C1,D1);
								baricentro.setLocation((solapamiento * (bar1.getX() + bar2.getX()))/solapamiento, (solapamiento * (bar1.getY() + bar2.getY()))/solapamiento);


							} else if(Y.getX() < xC1) {
								if(Q.getX() > xD1) {
									//*** Caso 18 ***
									//Calcular el area de NVR + NRY + VD1R
									double areaNVR = calculateArea(N,V,R);
									double areaNRY = calculateArea(N,R,Y);
									double areaVD1R = calculateArea(V,D1,R);
									solapamiento = solapamiento + areaNVR + areaNRY + areaVD1R;
									//Calcular el baricentro
									Point2D.Double bar1 = getBar(N,V,R);
									Point2D.Double bar2 = getBar(N,R,Y);
									Point2D.Double bar3 = getBar(V,D1,R);
									baricentro.setLocation((solapamiento * (bar1.getX()+bar2.getX()+bar3.getX()))/solapamiento, (solapamiento * (bar1.getY()+bar2.getY()+bar3.getY()))/solapamiento);

								} else if(Q.getX() < xD1) {
									if(xA2 < xB1) {
										//*** Caso 19 ***
										//Calcular el area de NVY + VYQ
										double areaNVY = calculateArea(N,V,Y);
										double areaVYQ = calculateArea(V,Y,Q);
										solapamiento = solapamiento + areaNVY + areaVYQ;
										//Calcular el baricentro 
										Point2D.Double bar1 = getBar(N,V,Y);
										Point2D.Double bar2 = getBar(V,Y,Q);
										baricentro.setLocation((solapamiento * (bar1.getX() + bar2.getX()))/solapamiento, (solapamiento * (bar1.getY() + bar2.getY()))/solapamiento);

									} else if(xA2 > xB1) {
										if(V.getX() > xA1) {
											//*** Caso 20 ***
											//Calcular WB1Y + WVY + VYQ
											double areaWB1Y = calculateArea(W,B1,Y);
											double areaWVY = calculateArea(W,V,Y);
											double areaVYQ = calculateArea(Y,V,Q);
											solapamiento = solapamiento + areaWB1Y + areaWVY + areaVYQ;
											//Calcular el baricentro
											Point2D.Double bar1 = getBar(W,B1,Y);
											Point2D.Double bar2 = getBar(W,V,Y);
											Point2D.Double bar3 = getBar(Y,V,Q);
											baricentro.setLocation((solapamiento * (bar1.getX()+bar2.getX()+bar3.getX()))/solapamiento, (solapamiento * (bar1.getY()+bar2.getY()+bar3.getY()))/solapamiento);

										} else if(V.getX() < xA1) {
											//*** Caso 21 ***
											//Calcular el area de A1YB1 + A1YQ
											double areaA1YB1 = calculateArea(A1,Y,B1);
											double areaA1YQ = calculateArea(A1,Y,Q);
											solapamiento = solapamiento + areaA1YB1 + areaA1YQ;
											//Calcular el baricentro
											Point2D.Double bar1 = getBar(A1,Y,B1);
											Point2D.Double bar2 = getBar(A1,Y,Q);
											baricentro.setLocation((solapamiento * (bar1.getX() + bar2.getX()))/solapamiento, (solapamiento * (bar1.getY() + bar2.getY()))/solapamiento);

										}
									}
								}
							}
						} else if(yA2 < yA1 && yB2 < yB1 && Y.getX() > xC1 && V.getX() > xD1) {
							//*** Caso 22 ***
							//Calcular el area de XA2Z + XC1Z
							double areaXA2Z = calculateArea(X,A2,Z);
							double areaXC1Z = calculateArea(X,C1,Z);
							solapamiento = solapamiento + areaXA2Z + areaXC1Z;
							//Calcular el baricentro
							Point2D.Double bar1 = getBar(X,A2,Z);
							Point2D.Double bar2 = getBar(X,C1,Z);
							baricentro.setLocation((solapamiento * (bar1.getX() + bar2.getX()))/solapamiento, (solapamiento * (bar1.getY() + bar2.getY()))/solapamiento);
						
						}


						double m1 = ((VEColi) a1).getMass();
						double m2 = ((VEColi) a2).getMass();
						double ba1 = ((VEColi) a1).getLength();
						double ha1 = ((VEColi) a1).getWidth();
						double lambdaG = 0.01;
						double lambdaD = 0.005; //0.025
						double k1 = (m2/m1)*lambdaD;
						double k2 = ( m2/m1 ) * (1/ (Math.pow(ba1, 2) + Math.pow(ha1, 2)))*lambdaG;

						double D = k1 * Math.sqrt(solapamiento);
						double G = k2 * solapamiento;

						double newHeading = 1;
						
						/**Estudio los casos en los que he hecho un movimiento simetrico anteriormente
						para ver como será el giro del angulo con el nuevo posicionamiento.**/
						if(cambio) {
							if(caso == 2) {
								if(baricentro.getX() * baricentro.getY() > 0 ){
									newHeading = ((VEColi) a1).getHeading() - G;
								}
							} else if(caso == 4) {
								if(baricentro.getX() * baricentro.getY() > 0 ){
									newHeading = ((VEColi) a1).getHeading() + G;
								} else
									newHeading = ((VEColi) a1).getHeading() - G;

							} else if(caso == 5) {
								if(baricentro.getX() * baricentro.getY() > 0){
									newHeading = ((VEColi) a1).getHeading() - G;
								} else 
									newHeading = ((VEColi) a1).getHeading() + G;
							} else if(caso == 6) {
								if(baricentro.getX() * baricentro.getY() > 0)
									newHeading = ((VEColi) a1).getHeading() - G;
							} else if(caso == 8) {
								if(baricentro.getX() * baricentro.getY() > 0)
									newHeading = ((VEColi) a1).getHeading() - G;
								else
									newHeading = ((VEColi) a1).getHeading() + G;
							}
						} else {
							if(caso == 1) {
								if(baricentro.getX() * baricentro.getY() > 0 ){
									newHeading = ((VEColi) a1).getHeading() - G;
								} else
									newHeading = ((VEColi) a1).getHeading() + G;

							} else if(caso == 3){
								if(baricentro.getX() * baricentro.getY() < 0 )
									newHeading = ((VEColi) a1).getHeading() + G;
							} else if(caso == 7) {
								if(baricentro.getX() * baricentro.getY() < 0 )
									newHeading = ((VEColi) a1).getHeading() - G;
							}

						}


						/*Movimiento lo consideramos descomponiendolo en dos partes
						 * 	1.Un desplazamiento de la bacteria en la direccion que une el centro de gravedad de las dos bacterias.
						 * 	2.Un giro que se determina segun la geometria del sistema.
						 */

						//Giro la bacteria.
						((VEColi) a1).setHeading( ((VEColi) a1).getHeading() + newHeading );
						//Lo movemos en la direccion del vector que va del centro de la bacteria 2 a la bacteria 1.
						//Desplazo la bacteria.
						
						double t = Math.sqrt(Math.pow(xG2, 2) + Math.pow(yG2,2));
						if(t != 0) {
							xG2 = xG2 /t;
							yG2 = yG2 /t;
						}
						
						space.moveTo(((VEColi) a1), space.getLocation(((VEColi) a1)).getX() + xG2*D, space.getLocation(((VEColi) a1)).getY() + yG2*D,0); 
						

					}
				}
			}
		}

	}//end Relaxation
} //end ColonyExpansion


/**
 * This class provides search colony edge  
 * TODO: Improve and make it a public class     
 */
class ColonyEdge {

	public static boolean isEdge(Object o, Iterable<Point2i> ring) {
		return (agentCount(o, ring) == 0);
	}

	public static int agentCount(Object o, Iterable<Point2i> ring) {
		//Context context = (Context) ContextUtils.getContext(o);
		//Grid grid = (Grid) context.getProjection("grid-space");
		int count= 0;

		for(Point2i p: ring) {
			count+= agentCount(o, p.getX(), p.getY());
		}
		return count;
	}

	@SuppressWarnings("rawtypes")
	public static int agentCount(Object o, int x, int y) {
		Context context = (Context) ContextUtils.getContext(o);
		Grid grid = (Grid) context.getProjection("grid-space");
		int count= 0;
		for (@SuppressWarnings("unused") Object a: grid.getObjectsAt(x, y)) {
			count++;
		}
		return count;
	}
}

/**
 * This class provides the iterable collection of 
 * colony ring coordinates
 * TODO: Improve and make it a public class     
 */
class ColonyRing {
	private int x= 0;
	private int y= 0;

	public ColonyRing() {
	}

	public ColonyRing(int xx, int yy) {
		x= xx;
		y= yy;
	}

	public Iterable<Point2i> Items(int r) {
		List<Point2i> list= new ArrayList<Point2i>();

		int rx= -r;
		int ry= 0;
		for(int i= 0; i < (r > 0 ? (r * 4) : 1); i++) {
			list.add(new Point2i(new int[] {x+rx,y+ry}));

			rx= (i % 2 == 0 ? rx + 1 : rx);
			ry= (i % 2 == 0 ? (r - Math.abs(rx)) : -(r - Math.abs(rx)));
		}
		Collections.shuffle(list);
		return list;
	}
}
