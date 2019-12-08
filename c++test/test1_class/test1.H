#include <iostream>
 
using namespace std;
 
class Line
{
   public:
      void setLength( double len );
      double getLength( void );
      double getDiameter( void );
      
       // constructors
	  Line(double len,double dia); 
	  // destructors 在每次删除所创建的对象时执行
	  //析构函数的名称与类的名称是完全相同的，
	  //只是在前面加了个波浪号（~）作为前缀，它不会返回任何值，也不能带有任何参数。
	  //析构函数有助于在跳出程序（比如关闭文件、释放内存等）前释放资源。
	  ~Line(); 
   private:
      double length;
      double diameter;
};
 
