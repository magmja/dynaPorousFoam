
vector exampleVector = vector::zero;

Info << "vector:" << tab << exampleVector << endl;

scalar listLength(5);

List<vector> exampleVectorList(listLength,vector::zero);

forAll(exampleVectorList, i)
{
    vector& v = exampleVectorList[i];

    v = vector(i,2*i,3*i);

    Info << "vector:" << tab << v << endl;
}

Info << endl << endl;

List<List<vector> > exampleVectorListList(listLength,exampleVectorList);

forAll(exampleVectorListList,i)
{
    List<vector>& vList = exampleVectorListList[i];

    forAll(vList,j)
    {
        vector& v = vList[j];

        Info << "vector:" << tab << v << endl;
    }
}