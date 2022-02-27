/**
 * 
 *     Rangement de Cuthill Mc Kee  
 *     et
 *     Resolution avec LDLT
 *     -------
 * 
 *     RATEFIARISON
 *     Harivony Lalatiana
 *     MISA - M1
 *      
 * 
 **/


#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>

using namespace std;

/***************    BEGIN  LDLT    *********************/

int dim(10);            // definir la taille du systeme 

/*
    Allocation/Affichage de matrice et de vecteur
*/

float **newMat(int rows, int cols);     // creer une matrice rows x cols
float *newVec(int n);                   // creer un vecteur de taille n
void displayMat(float** mat);           // afficher une matrice de mat
void displayVec (float* v);             // afficher un vecteur
void displayVec (int* v);             // afficher un vecteur

/*
 *
 *  CLASS: LDLt
 * 
 * 
*/

class LDLT{
// Methods
public:
    LDLT();
    ~LDLT();
    void loadAb();
    void computeP_i();
    void computeAp_nDiag();
    void computeAP();
    void computeX();

    // getter
    float** getA(){return A;}
    float* getb(){return b;}
    float** getD(){return D;}
    float* getX(){return x;}
    float getAp(int i, int j);
    void setAp(int i, int j, float val);

    void solveTriangInf();
    void solveDiag();
	void solveTriangSup();

    void profil(float** mat, string name);
    void getProfil(float* mat, string name);
    int* solve();
    void gaussSolve();
    vector<float> AP = {};
	vector<int> l = {};
	vector<int> p = {};
    vector<int> nDiag;

// Attributs
private:
    float** A;
    float** L;
    float** D;
    float* x;    // la solution

    float* b;    // second membre
    float* p_i;  // indice des premieres element non null
};


/***************    END  LDLT    *********************/

class Graphe {
    public: 
        Graphe(string filename);
        Graphe(vector<int>* neighbors, int dim);
        Graphe();
        ~Graphe();

        //les getters et setters
        int getDim(){return _dim;}
        vector<int>* getNeigborsList(){return _neighbors;}
        vector<int> getLastN();
        void setFilename(string filename);

        //gerer les affichages
        void displayNeighborsList();
        void displayVect1D(vector<int> vect);
        void displayVect2D(vector<vector<int>> vect);
        
        int findEccentricity(int node);//retourne l'excentricité d'un noeud

    
    private:
        // void findNeighbors(int node, vector<vector<int>> N); //retourne la liste des voisins d'un noeud

        string _file; //fichier qui contient la matrice pour former la graphe
        int _dim; //dimension de la matrice de la graphe et aussi le nombres de noeuds de la graphe
        vector<int> _lastNlist; //Contient la liste de N avant le N = {}
        vector<int>* _neighbors; //Contient la liste des voisinages de tous les noeuds
};

class CuthillMackee
{
    public:
        CuthillMackee(string filename);
        ~CuthillMackee();

        int** getP(){return _P;}

        void solve(int node);
        void storeData(string filename);
        //Les opérations matricielles
        int* matTimesVect(int** Mat, int* vect, int row1, int col1, int vectrow);//Produit d'une matrice à un vecteur
        int** matTimesMat(int** Mat1, int** Mat2, int row1, int col1, int row2, int col2);//Produit de deux matrices
        int** transpose(int** M, int row, int col); //Transposer d'une matrice
        void displayArray(int* arr, int row);
        void displayMatrix(int** Mat, int row, int col);

    private:
        bool is_inSigma(int node);
        bool isIn(vector<vector<int>> vect, int node);
        vector<int> sortVect(vector<int> vect);
        vector<vector<int>> minNeighborSortedList(vector<int>* neighbors, int* sigma, vector<int> node);
        void findFirstNode(int node); //Recherche de 1er noeud
        void buildSigma();
        void buidSigmaInverse();
        void getData();
        void buildP(int* sigma); //Chercher la matrice de passage P;
        
        

        string _filename;
        Graphe* _graphe1;
        int _firstNode;
        int _smNode;//Noeud quelconque qu'on choisit, on le stocke pour remplir sigma
        int _dim;
        int* _sigma;
        int* _b;
        int** _P;
        int** _A;

};

int main(){
    
    CuthillMackee* cm = new CuthillMackee("matrixSystem.txt");
    cm->solve(3);//Solve prend un sommet (int)
    cm->storeData("matWithMinProfile.txt");

    LDLT ldlt;// factorisation et resolution
    
    cout << "____________________________________" << endl;
    cout << "\n  Calcule de la solution final    " << endl;
    cout << "____________________________________" << endl << endl;
    cout << "On a A'.x' = b' ,  \n => Pt.x' = b' \n => x = P.x'"<< endl<< endl;
    // matTimesVect(transpose(_P, _dim, _dim), _b, _dim, _dim, _dim); //Pt*b
    
    
    int* solution = new int [10];
    solution = ldlt.solve();
    
    cout << "____________________________________" << endl << endl;
    cout << "       La solution x finale est      "<< endl;
    cout << "____________________________________" << endl << endl;
    
    // x = P.x'
    solution = cm->matTimesVect(cm->getP(),solution,dim,dim,dim);
    
    displayVec(solution);
    
    
    // succes de l'execution
    return 0;
}


/**
 * sans parametre
 ******************/
Graphe::Graphe(){}

/**
 * param: nom de fichier où se trouve 
 * la dimension et la matrice dont on formera la graphe
 * *****************************************************/
Graphe::Graphe(string filename)
{
    _file = filename;
    ifstream file(_file);

    if (file.is_open())
    {
        cout << "Ouverture de " << _file << " pour recuperer la liste des voisinages :" << endl;
        
        file >> _dim;// dimension de la matrice
        string temp("");
        _neighbors = new vector<int>[_dim];

        for (size_t j = 0; j < _dim; j++)
        {
            for (size_t i = 0; i < _dim; i++)
            {
                getline(file, temp, ',');
                if ( (atof(temp.c_str()) != 0) && (i != j) )
                {
                    _neighbors[j].push_back(i);
                }
            }
        }

    }
}

/**
 * param: tableau de vecteur
 * param: dimension de la matrice dont on veut contruire sa graphe
 *****************************************************************/
Graphe::Graphe(vector<int>* neighbors, int dim)
{
    _file = "";
    _dim = dim;
    _neighbors = new vector<int>[_dim];

    for (size_t j = 0; j < _dim; j++)
    { 
        _neighbors[j] = neighbors[j];
    }
}


/**
 * Quand on change le nom de fichier, cela signifie qu'on
 * a changer de matrice. Donc, on change les dimensions, 
 * les noeuds de la graphe. (c-à-d: _file, _dim, _neighbors)
 ***************************************/
void Graphe::setFilename(string filename)
{
    _file = filename;
    ifstream file(_file);

    if (file.is_open())
    {
        cout << "File opening with success :-)" << endl;
        
        file >> _dim;// dimension de la matrice
        string temp("");
        _neighbors = new vector<int>[_dim];

        for (size_t j = 0; j < _dim; j++)
        {
            for (size_t i = 0; i < _dim; i++)
            {
                getline(file, temp, ',');
                if ( (atof(temp.c_str()) != 0) && (i != j) )
                {
                    _neighbors[j].push_back(i);
                }
            }
        }
    }
}


Graphe::~Graphe()
{
    delete[] _neighbors;
}

void Graphe::displayNeighborsList()
{
    cout << "==============================================" << endl;
    cout << "Liste des voisinages des noeuds de la graphe :" << endl;
    cout << "==============================================" << endl;
    for (size_t i = 0; i < _dim; i++)
    {
        cout << "[" << i << "] : ";
        for (size_t j = 0; j < _neighbors[i].size(); j++)
        {
            cout << _neighbors[i][j] << "  ";
        }
        cout << endl;
    }
}

/**
 * param: int, numero noeud
 * return: int, excetricité du noeud
 ****************************************************/
int Graphe::findEccentricity(int node)
{
    vector<int> Ncurrent;
    vector<int> Nprev;
    vector<int> tempv;
    int it(1); //profondeur du niveau
    int itemp(0); //index temporaire

    cout <<"==============" << endl;
    cout << "    n = " << node << endl;
    cout <<"==============" << endl;

    //N0 et N1
    Nprev.push_back(node);
    for (size_t i = 0; i < _neighbors[node].size(); i++)
    {
        Ncurrent.push_back(_neighbors[node][i]);
    }
    cout << "N[0] = {" ; displayVect1D(Nprev); cout << "}" << endl;    
    
    while (Ncurrent.size() != 0)
    {
        cout << "N[" << it << "] = {" ; displayVect1D(Ncurrent); cout << "}" << endl;
        it++;

        //Chercher les voisins de Ncurrent
        for (size_t i = 0; i < Ncurrent.size(); i++)
        {
            itemp = Ncurrent[i];
            for (size_t j = 0; j < _neighbors[itemp].size(); j++)
            {
                tempv.push_back(_neighbors[itemp][j]);
            } 
        }
        
        //Après avoir eu voisin des élements de Ncurrent
        //1)Eliminons les doublons dans tempv
        for (size_t i = 0; i < tempv.size(); i++)
        {
            for (size_t j = 0; j < i; j++)
            {
                if (tempv[i] == tempv[j])
                {
                    tempv.erase(tempv.begin() + i);
                    i--;
                }
            } 
        }
        
        //2)Eliminons les éléments de Ncurrent dans tempv
        for (size_t i = 0; i < tempv.size(); i++)
        {
            for (size_t j = 0; j < Ncurrent.size(); j++)
            {
                if (tempv[i] == Ncurrent[j])
                {
                    tempv.erase(tempv.begin() + i);
                    i--;
                }
            } 
        }
       
       //3)Eliminons les éléments de Nprev dans tempv
       for (size_t i = 0; i < tempv.size(); i++)
        {
            for (size_t j = 0; j < Nprev.size(); j++)
            {
                if (tempv[i] == Nprev[j])
                {
                    tempv.erase(tempv.begin() + i);
                    i--;
                }
            } 
        }

        Nprev = Ncurrent;
        Ncurrent = tempv;
        tempv.clear();
    }
    
    if (Ncurrent.size() == 0)
    {
        _lastNlist = Nprev;
        cout << "N[" << it << "] = {}" << endl;
    }
    

    return (it - 1);
}


void Graphe::displayVect1D(vector<int> vect)
{
    for (int i = 0; i < vect.size(); i++)
    {
        cout << vect[i] << "   ";
    } 
}


void Graphe::displayVect2D(vector<vector<int>> vect)
{
    for (int i = 0; i < vect.size(); i++)
    {
        for (int j = 0; j < vect[i].size(); j++)
        {
            cout << vect[i][j] << "  ";
        }
        cout << endl;
    } 
}

vector<int> Graphe::getLastN()
{
    return _lastNlist;
}



/**
 * param: string, nom de fichier
 *******************************/
CuthillMackee::CuthillMackee(string filename)
{
    _filename = filename;
    _graphe1 = new Graphe(_filename);
    getData();
}

/**
 * Récuperer les valeurs de A et b
 * dans un fichier appellée filename
 **********************************/
void CuthillMackee::getData()
{
    ifstream file(_filename);
    if (file.is_open())
    {
        cout << "Ouverture de " << _filename << " pour récuperer la matrice A et le vecteur b :" << endl;
        string temp("");

        file >> _dim;
        _b = new int[_dim];
        _A = new int*[_dim];
        _P = new int*[_dim];
        for (size_t i = 0; i < _dim; i++)
        {
            _A[i] = new int[_dim];
            _P[i] = new int[_dim];
        }
        _sigma = new int[_dim];


        for (size_t i = 0; i <= _dim; i++)
        {
            for (size_t j = 0; j < _dim; j++)
            {
                if (i == _dim)
                {
                    getline(file, temp, ',');
                    _b[j] = atof(temp.c_str());
                    // cout << _b[j] << "  ";
                }
                else
                {
                    getline(file, temp, ',');
                    _A[i][j] = atof(temp.c_str());
                    // cout << _A[i][j] << "  ";
                }
            }
            // cout << endl;
        }
    }
}


void CuthillMackee::solve(int node)
{
    


    /**
     * Etape de cuthill-Mackee
     * *********************************/
    cout <<"\n__________________________________________________________" << endl;
    cout << " Voici la matrice matrice originale avec son second membre :" << endl;
    cout << "____________________________________________________________" << endl;
    displayMatrix(_A, _dim, _dim);
    cout << endl;
    displayArray(_b, _dim);
    cout << endl;
    
    
    findFirstNode(node);
    buildSigma();
    buildP(_sigma);//sigma

    cout <<"\n-------------------------------------------------------------" << endl;
    cout << "           La matrice de passage trouver                   :" << endl;
    cout <<"-------------------------------------------------------------" << endl;
    displayMatrix(_P, _dim, _dim);
    cout << endl;

    
    //Calcul de Aprim et de bprim
    _A = matTimesMat(transpose(_P, _dim, _dim), _A, _dim, _dim, _dim, _dim); //Pt*A
    _A = matTimesMat(_A, _P, _dim, _dim, _dim, _dim);//A*P
    _b = matTimesVect(transpose(_P, _dim, _dim), _b, _dim, _dim, _dim); //Pt*b

    cout <<"------------------------------------------------------------------------------" << endl;
    cout << " (cuthill-Mackee) Voici la matrice matrice optimisée avec son second membre :" << endl;
    cout <<"------------------------------------------------------------------------------" << endl;


    displayMatrix(_A, _dim, _dim);
    cout << endl;
    displayArray(_b, _dim);
    cout << endl;

    /**
     * Etape de cuthill-Mackee inverse
     * *********************************/
    buidSigmaInverse();
/*
    //Calcul de Aprim et de bprim
    _A = matTimesMat(_P, _A, _dim, _dim, _dim, _dim); //P*A
    _A = matTimesMat(_A, transpose(_P, _dim, _dim), _dim, _dim, _dim, _dim);//A*Pt
    _b = matTimesVect(_P, _b, _dim, _dim, _dim); //P*b
    
    cout << "Debug P" << endl;
    displayMatrix(_P, _dim, _dim);

    cout <<"--------------------------------------------------------------------------------------" << endl;
    cout << " (cuthill-Mackee-Inverse) Voici la matrice matrice optimisée avec son second membre :" << endl;
    cout <<"--------------------------------------------------------------------------------------" << endl;
    displayMatrix(_A, _dim, _dim);
    cout << endl;
    displayArray(_b, _dim);
    cout << endl;
    */
}

/*
*   Recherche du premier node
*/
void CuthillMackee::findFirstNode(int node)
{
    cout << "Choix du premier noeud :" << endl;
    int excetricity;
    _firstNode = node;
    _smNode = node;
    vector<int> vect;
    
    excetricity = _graphe1->findEccentricity(node);
    vect = _graphe1->getLastN();
    
    for (size_t i = 0; i < vect.size(); i++)
    {
        if (excetricity < _graphe1->findEccentricity(vect[i]))
        {
            _firstNode = vect[i];
            excetricity = _graphe1->findEccentricity(vect[i]);
        }
    }

    cout <<"Premier noeud trouver est : " << _firstNode << endl;
}

/*
*   Construction de sigma
*/

void CuthillMackee::buildSigma()
{
    cout <<"---------------------" << endl;
    cout << "Recherche de sigma : " << endl;
    cout <<"---------------------" << endl;
    int it(0);
    vector<int> vect = {_firstNode}; //lister dans cet vecteur, par ordre croissant, les noeuds le plus petit voisin non numéroté
    vector<vector<int>> vectTemp; //vecteur temporaire de vect


    //Initialiser sigma pour vérifier les noeuds déja numérotés
    _sigma[0] = _firstNode;
    for (size_t i = 1; i < _dim; i++)
    {
        _sigma[i] = -1;
    }

    //remplir sigma
    while(vect.size() != 0)
    {
        vectTemp = minNeighborSortedList(_graphe1->getNeigborsList(), _sigma, vect);
        vect.clear();
        for (size_t i = 0; i < vectTemp.size(); i++)
        {
            for (size_t j = 0; j < vectTemp[i].size(); j++)
            {
                vect.push_back(vectTemp[i][j]);
            }
        }

        for (size_t i = 0; i < vect.size(); i++)
        {
            for (size_t j = 0; j < _dim; j++)
            {
                if (_sigma[j] == vect[i])
                {
                    break;
                }
                else if (_sigma[j] == -1)
                {
                    _sigma[j] = vect[i];
                    break;
                }
            }
        }
    }

    //affichage sigma
    for (size_t i = 0; i < _dim; i++)
    {
        cout << _sigma[i] << "   ";
    }
    cout << endl;
}

/**
 * sigmaInverse[i] = (dim-1) - sigma[i]
 * dim - 1 car on part de l'indice 0
 **************************************/
void CuthillMackee::buidSigmaInverse()
{
    cout <<"-----------------------------" << endl;
    cout << "Recherche de sigmaInverse : " << endl;
    cout <<"-----------------------------" << endl;

    for (size_t i = 0; i < _dim; i++)
    {
        _sigma[i] = (_dim-1) - _sigma[i];
    }

    for (size_t i = 0; i < _dim; i++)
    {
        cout << _sigma[i] << "   ";
    }
    cout << endl;
}


/**
 * params: .vector<int>*, liste des voisins des noeuds de la graphe(c'est la graphe)
 *         .vector<int>, liste des éléments de sigma
 *         .vector<int>, liste de numéros des noeud(id des noeud)
 * 
 * return: liste des voisins non numérotés des noeuds de node. (rangé dans l'ordre croissant
 *         suivant le nombre de voisins non numérotés des noeud de node)
 * 
 ************************************************************************************************/
vector<vector<int>> CuthillMackee::minNeighborSortedList(vector<int>* neighbors,int* sigma, vector<int> nodes)
{
    int iTemp;
    int currentNode;
    int counter(0); //compte le nombre de voisins non numerotés d'un noeud(des voisins de node)
    // vector<int> vectTemp;
    vector<vector<int>> notNumeroted;
    vector<vector<int>> vCounter; //stocke les counter
    vector<vector<int>> vCounterSorted; //stocke les counter ordonnées dans l'ordre croissant
    
    for (size_t it = 0; it < nodes.size(); it++)
    {
        vector<int> vTemp;
        notNumeroted.push_back(vTemp);
        vCounter.push_back(vTemp);
        vCounterSorted.push_back(vTemp);

        //Chercher les voisins non numérotés de node[it]
        currentNode = nodes[it];
        for (size_t i = 0; i < neighbors[currentNode].size(); i++)
        {
            // cout << neighbors[currentNode][i] << endl;
            if ( !is_inSigma(neighbors[currentNode][i]) )
            {
                // vTemp.push_back(neighbors[currentNode][i]);
                notNumeroted[it].push_back(neighbors[currentNode][i]);
                for (size_t k = 0; k < neighbors[neighbors[currentNode][i]].size(); k++)
                {
                    if (!is_inSigma(neighbors[neighbors[currentNode][i]][k]) && !isIn(notNumeroted, neighbors[neighbors[currentNode][i]][k]))
                    {
                        counter++;
                    }
                }
                // cout <<"Counter : " << counter << endl;
                vCounter[it].push_back(counter);
                counter = 0;
            }
            
        }
        // for (size_t l = 0; l < vTemp.size(); l++)
        // {
        //     notNumeroted[it].push_back(vTemp[l]);
        // }
        // vTemp.clear(); 

        vCounterSorted[it] = sortVect(vCounter[it]);

        //Ranger les voisins non numérotés de node[it] par rapport à leurs nombres de ces voisins non numérotés
        for (size_t i = 0; i < vCounterSorted[it].size(); i++)
        {
            for (size_t j = 0; j < vCounter[it].size(); j++)
            {
                //on peut retrouver notNumeroted à l'aide de vCounter[it] car ils ont le meme indice
                //c-à-d que vCounter[it][1] est le nombre de voisins non numéroté du noeud notNumeroted[1]
                if (vCounter[it][j] == vCounterSorted[it][i])
                {
                    iTemp = notNumeroted[it][i];
                    notNumeroted[it][i] = notNumeroted[it][j];
                    notNumeroted[it][j] = iTemp;

                    vCounter[it][i] = -1;//Pour éviter de repeter la permutation sur ce i èm élement
                    break;
                }
            }
        }
    }

    return notNumeroted;
}


void CuthillMackee::buildP(int* sigma)
{
    for (size_t i = 0; i < _dim; i++)
    {
        for (size_t j = 0; j < _dim; j++)
        {
            if (i == sigma[j])
            {
                _P[i][j] = 1;
            }
            else
            {
                _P[i][j] = 0;
            }
        }
    }
}



void CuthillMackee::storeData(string filename)
{
    ofstream outFile(filename);

    outFile << to_string(_dim) << endl;

    for (size_t i = 0; i <= _dim; i++)
    {
        for (size_t j = 0; j < _dim; j++)
        {
            if (i==_dim)
            {
                outFile << to_string(_b[j]) << " ";
            }
            else
            {
                outFile << to_string(_A[i][j]) << " ";
            }
        }
        outFile << endl;
    }  
}


/**
 * Une fonction qui vérifie l'existance d'un noeud dans sigma
 */
bool CuthillMackee::is_inSigma(int node)
{


    for (size_t i = 0; i < _dim; i++)
    {
        if (node == _sigma[i])
        {
            return true;
        }
    }
    
    return false;
}

/**
 * C'est une fonction pour verifier l'appartenance
 * d'un int dans l'éléments d'un vecteur de vecteur
 ************************************************/
bool CuthillMackee::isIn(vector<vector<int>> vect, int node)
{
    for (size_t i = 0; i < vect.size(); i++)
    {
        for (size_t j = 0; j <vect[i].size(); j++)
        {
            if (node == vect[i][j])
            {
                return true;
            }
        }
    }
    
    return false;
}


/**
 * Arranger dans l'ordre croissant un vecteur
 * param: vector<int>, vecteur
 * return: vector<int>, vecteur
 *******************************/
vector<int> CuthillMackee::sortVect(vector<int> vect)
{
    int temp;
    for (size_t i = 0; i < vect.size(); i++)
    {
        for (size_t j = 0; j < vect.size(); j++)
        {
            if (vect[j] < vect[j-1])
            {
                temp = vect[j-1];
                vect[j-1] = vect[j];
                vect[j] = temp;
            }
        }
        
    }
    
    return vect;
}

/**
 * Produit de deux matrices
 * param: int**, matrice 1
 *        int**, matrice 2
 *        dimensions des deux matrices (dim de la mat1 et puis mat2)
 * return: int **
 *******************************/
int** CuthillMackee::matTimesMat(int** Mat1, int** Mat2, int row1, int col1, int row2, int col2)
{
    //T va contenir le produit des deux matrices
    //T a comme dimension (row1, col2)
    int** T = new int*[row1];
    for (size_t i = 0; i < row1; i++)
    {
        T[i] = new int[col2];
        for (size_t j = 0; j < col2; j++)
        {
            T[i][j] = 0;
        }
        
    }


    if (col1 == row2)
    {   
        for (size_t k = 0; k < col2; k++)
        {
            for (size_t i = 0; i < row1; i++)
            {
                for (size_t j = 0; j < col1; j++)
                {
                    T[i][k] += Mat1[i][j] * Mat2[j][k];
                }  
            }
        } 
    } 
    else
    {
        cout << "Erreur: Nombre de col = " << col1 << " est different nombre de ligne = " << row2 << " !" << endl;
    }

    
    return T;
}

/**
 * Produit d'une matrice avec un vecteur
 * param: int**, matrice
 *        int*, vecteur
 *        int, dimension de la matrice (row1, col1)
 *        int, dimension du vecteur (vectrow)
 * return: int*, vecteur
 *****************************************/
int* CuthillMackee::matTimesVect(int** Mat, int* vect, int row1, int col1, int vectrow)
{
    int* T = new int[vectrow];
    for (size_t i = 0; i < vectrow; i++)
    {
        T[i] = 0;
    }


    if (col1 == vectrow)
    {   
        for (size_t i = 0; i < row1; i++)
        {
            for (size_t j = 0; j < col1; j++)
            {
                T[i] += Mat[i][j] * vect[j];
            }  
        }
    } 
    else
    {
        cout << "Erreur: Nombre de col = " << col1 << " est different nombre de ligne = " << vectrow << " !" << endl;
    }

    return T;
}

/**
 * Transposer d'une matrice
 * param: int**, la matrice à transposer
 *        int, la dimension de matrice (line, colonne)
 * 
 * return: int**, la matrice transposée
 ***********************************************************/
int** CuthillMackee::transpose(int ** Mat, int row, int col)
{
    int** T = new int*[row];
    for (size_t i = 0; i < row; i++)
    {
        T[i] = new int[col];
        for (size_t j = 0; j < col; j++)
        {
            T[i][j] = 0;
        }
        
    }

    for (size_t i = 0; i < row; i++)
    {
        for (size_t j = 0; j < col; j++)
        {
            T[i][j] = Mat[j][i];
        }
        
    }
    
    return T;
}

void CuthillMackee::displayArray(int* arr, int row)
{
    for (size_t i = 0; i < row; i++)
    {
        cout << arr[i] << "  ";
    }
    cout << endl;
}

void CuthillMackee::displayMatrix(int** Mat, int row, int col)
{
    for (size_t i = 0; i < row; i++)
    {
        for (size_t j = 0; j < col; j++)
        {
            cout << Mat[i][j] << "  ";
        }
        cout << endl;        
    }
    
}


/***************    BEGIN  LDLT    *********************/

/***** MATRICE - VECTEUR ******/

float **newMat(int rows, int cols){ //allocation dynamique de tableau rows*cols
    float **mat(NULL);
    mat = new float*[rows];

    for (int i=0;i<rows;i++)        
        mat[i] = new float[cols];
    return mat;
}

float *newVec(int n){               //alllouer un vecteur de dim n
    float* b;
    b  = new float [n];
    return b;
}

void displayMat(float** mat){    //afficher la matrice
    for (int i=0;i<dim;i++){
        for(int j=0; j<dim; j++)
            cout << mat[i][j] << " ";
        cout << endl;
    }
}

void displayVec (float* v){         //afficher un vecteur
    for(int i=0; i<dim; i++){
        cout << "[" << v[i] << "]" << endl;
    }
}

void displayVec (int* v){         //afficher un vecteur
    for(int i=0; i<dim; i++){
        cout << "[" << v[i] << "]" << endl;
    }
}


/*
 *
 *  DEFINITION : METHODE LDLT
 * 
 * 
*/

LDLT::LDLT(){
    // Allcation des matrice et des vecteurs
    A = newMat(dim,dim);
    L = newMat(dim,dim);
    D = newMat(dim,dim);
    b = newVec(dim);
    x = newVec(dim);
}

// Chargger la matrice depuis un fichier
void LDLT::loadAb(){
    string fileName;
    ifstream inFile;
    
    //cout << "Entrez le fichier a charger : "; cin >> fileName;

    // ouverture du fichier, 
    inFile.open("matWithMinProfile.txt");

    if(inFile){
        //recuperer la taille de la matrice
        inFile >> dim;

        // recuperer la partie inférieure avec le duagonal de la matrice
        for(int i=0; i<dim;i++){
            for(int j=0; j<dim;j++){
                inFile >> A[i][j];
            }
        }  
        
        // recuperer le second membre
        for(int i=0; i<dim;i++){
                inFile >> b[i];
        }
    }
    
    // Alerter l'erreur d'ouverture du fichier
    else
        cout << "erreur de lecture" << endl;
}

// charger le profil p_i de A
void LDLT::computeP_i(){
    p_i = newVec(dim);

    for(int i=0;i<dim;i++){
        bool logicNull = true;
        for(int j=0;j<=i; j++){
            if(A[i][j]==0){
                if(logicNull)
                    cout <<"  ";
                else
                    cout <<". ";
            }
            else {
                cout << "x ";
                if(logicNull)
                    p_i[i] = j;        // enregistrer la valeur de j comme etant le premmier element non null de la ligne
                logicNull = false;     // les zeros qui suive ne sont plus null logiquement
            }
        }
        cout << endl;
    }
    cout << endl;
}


// afficher le profi de la matrice mat
void LDLT::profil(float** mat, string name){
    cout << "Profil de la matrice " << name << " :"<< endl;
    for(int i=0;i<dim;i++){
        // eviter de faire un boucle de plus
        bool logicNull = true;
        for(int j=0;j<=i; j++){
            if(mat[i][j]==0){
                if(logicNull)
                    continue;
                else
                    cout <<"("<< i+1 << ";" << j+1 << "), ";
            }
            else {
                cout <<"("<< i+1 << ";" << j+1 << "), ";
                logicNull = false; // les zeros qui suive ne sont plus null logiquement
            }
        }
        cout << endl;
    }
    cout << endl;
}

void LDLT::computeAp_nDiag(){
	ifstream file("matWithMinProfile.txt", ios::in);
    if (file){
        file >> dim;

        // Initialisations depuis le fichier
        int ji(0);
        int count(-1); // Une variable qui va numeroter les éléments du profil
        float elt(0);
        char car('c');

        for(int_fast16_t i=0; i<dim; i++){
            ji = 0; 
            elt = 0;

            // Prendre seulement les éléments du profil
            for(int_fast16_t j=0; j<=i; j++){
                file >> elt;
                if (elt != 0 || ji != 0){
                    ji++;
                    if(ji == 1){
                        p.push_back(j);
                        l.push_back(i-j);
                    }
                }
                if(ji){
                    count++;
                    AP.push_back(elt);
                    if( i == j){
                        nDiag.push_back(count);
                    }
                }
            }            
            // Boucler jusqu'au dernier élément de la ligne
            while(car != '\n'){
                file.get(car);
            }
            car = 'c';
        }

        for(int_fast16_t i=0; i<dim; i++){
            file >> b[i];
        }
    }
}

// determiner la valeur de L et D avec A = L.D.Lt
void LDLT::computeAP(){     
    float sum = 0;
    for (int i=0; i<dim; i++){
        for(int j=0; j<i; j++){
            sum = 0;
            for(int k=0; k<j; k++){
                if (k>=p[i] && k>= p[j])
                    sum += AP[nDiag[i]-i+k]*AP[nDiag[k]]*AP[nDiag[j]-j+k];
            }
            if(j>=p[i] && j>=p[j])
                AP[nDiag[i]-i+j] = (1/AP[nDiag[j]])*(AP[nDiag[i]-i+j] - sum);
        }
        sum = 0; 
        for (int k=0; k<i; k++){
            if (k>=p[i])
                sum += AP[nDiag[k]]*(AP[nDiag[i]-i+k]*AP[nDiag[i]-i+k]);
        }
        AP[nDiag[i]] = AP[nDiag[i]] - sum;
    }
}

void LDLT::computeX(){
    // D'abord Lx = b
    float sum(0);
    for(int i=0; i<dim; i++){
        sum = 0;
        for(int j=0; j<i; j++){
            if(j>=p[i])
                sum += AP[nDiag[i]-i+j]*x[j];
        }
        x[i] = b[i] - sum;
    }
    // Puis Dx = x
    for(int i=0; i<dim; i++){
        x[i] = (1/AP[nDiag[i]])*x[i];
    }
    // Enfin Ltx = x
    for(int i= int(dim-1); i>=0; i--){
        sum = 0;
        for(int j=i+1; j<dim; j++){
            if(i>=int(p[j]))
                sum += AP[nDiag[j]-j+i]*x[j];
        }
        x[i] = x[i]-sum;
    }
}

float LDLT::getAp(int i, int j){
	
	if(j >= p[i] && j<=i){
		return AP[nDiag[i] - i + j];
	}	
	else 
	return 0;
}

void LDLT::setAp(int i, int j, float val){
	if(j >= p[i]) AP[nDiag[i] - i + j] = val;
}

// Appeler par etape les methodes de resolutions
int* LDLT::solve(){
    // 1 - charger les donner dempuis les fichier
    loadAb();
    
    // Interface
    cout << "La matrice A' :" << endl;
    displayMat(A);
    
    cout << "\nLe second membre b' :" << endl;
    displayVec(b);

    cout << "\n'Profil de A' :" << endl;
    computeP_i(); // tracer le profil

    // 2 - calcule de nDiag, l et p
    computeAp_nDiag();
    // 3 - Ajouter les Valeur de AP 
    computeAP(); 
    
    // 4 - Calculen de la solution
    computeX();  

    cout << "La solution (x') du systeme est :"<< endl;
    displayVec(x);                   // afficher le resultat
    
    int* xVect = new int[dim];
    for(int i=0;i<dim;i++)
        xVect[i] = (float) x[i];

    return xVect;                   // renvyer le resultat
}

// liberation de memoire apres execution
LDLT::~LDLT(){
	delete[] x;
	delete[] b;
    l.clear();
    p.clear();

    AP.clear();
    nDiag.clear();
	for(int i=0; i< dim; i++){
        delete[] D[i];
        delete[] A[i];
        delete[] L[i];
    }
    delete[] D;
    delete[] A;
    delete[] L;
}

/***************    END  LDLT    *********************/