#include "test1.H"
// constructors
Line::Line(double len, double dia)
{	
    cout << "Object is being created, length is not defined = " << length << endl;
    length = len;
    diameter = dia+0.2;
    cout << "Object is being created, length is defined= " << length << endl;
}
Line::~Line()
{
    cout << "Object is being deleted" << endl;
}
// member function
void Line::setLength( double len )
{
    length = len;
}
 
double Line::getLength( void )
{
    return length;
}

double Line::getDiameter( void )
{
    return diameter;
}


// main
int main( )
{
   Line line(10.0,0.1);
	cout << "Length of line : " << line.getLength() <<endl;

   line.setLength(6.0); 
   cout << "Length of line : " << line.getLength() <<endl;
 
   cout << "Diameter of line : " << line.getDiameter() <<endl;
   if (3)
   {Line line2(50.0,5.1);
	cout << "Length of line2 : " << line2.getLength() <<endl;
    cout << "Diameter of line2 : " << line2.getDiameter() <<endl;
   }
	   cout << "Length of line : " << line.getLength() <<endl;
 
   cout << "Diameter of line : " << line.getDiameter() <<endl;
   return 0;
}
