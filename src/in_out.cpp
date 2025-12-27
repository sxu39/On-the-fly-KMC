#include "in_out.h"
#include <iomanip>
using namespace KMC;
using namespace std;

namespace KMC
{
    static vector<prec_t> operator*(const vector<prec_t> &crd, const array<prec_t, 9> &box)
    {
        try
        {
            if (crd.size() != 3)
                throw std::invalid_argument("The size of a coordinate must be 3.");
            vector<prec_t> res(crd.size(), 0);
            for (int i = 0 ; i < 3 ; ++i)
                for (int j = 0 ; j < 3; ++j)
                    res[i] += crd[j] * box[3*j+i];
            return res;
        }
        catch(const std::exception& e)
        {
            cerr << e.what() << "\n";
        }
    }

    static prec_t det(const array<prec_t, 9> &mat)
    {
        return	mat[0]*mat[4]*mat[8] - mat[0]*mat[7]*mat[5] + mat[3]*mat[7]*mat[2] - 
                mat[3]*mat[1]*mat[8] + mat[6]*mat[1]*mat[5] - mat[6]*mat[4]*mat[2];
    }

    static array<prec_t, 9> inverse(const array<prec_t, 9> &mat)
    {
        prec_t d = det(mat);
        return array<prec_t, 9>{(mat[4]*mat[8] - mat[5]*mat[7]) / d, -(mat[1]*mat[8] - mat[2]*mat[7]) / d,
                                (mat[1]*mat[5] - mat[2]*mat[4]) / d, -(mat[3]*mat[8] - mat[5]*mat[6]) / d,
                                (mat[0]*mat[8] - mat[2]*mat[6]) / d, -(mat[0]*mat[5] - mat[2]*mat[3]) / d,
                                (mat[3]*mat[7] - mat[4]*mat[6]) / d, -(mat[0]*mat[7] - mat[1]*mat[6]) / d,
                                (mat[0]*mat[4] - mat[1]*mat[3]) / d};
    }
}

KMC::in_out::in_out():out_log("LOG"), out_phy("ENERGY"), quan_wid(16)
{
    out_phy.setf(ios::fixed);
    out_phy.setf(ios::showpoint);
    out_phy.setf(ios::left);
}

KMC::in_out::in_out(const std::string &log_file, const std::string &phy_file, size_t quan_wid):out_log(log_file), out_phy(phy_file), 
quan_wid(quan_wid)
{
    out_phy.setf(ios::fixed);
    out_phy.setf(ios::showpoint);
    out_phy.setf(ios::left);
}

KMC::in_out::~in_out()
{
    out_log.close();
    out_phy.close();
}

void KMC::in_out::set_log(const std::string &log_file)
{
    out_log.close();
    out_log.open(log_file);
}

void KMC::in_out::set_phy(const std::string &phy_file)
{
    out_phy.close();
    out_phy.open(phy_file);
}

const cell KMC::in_out::read_stru(const std::string &stru_file)
{
    ifstream in_stru(stru_file, ios::in);
    array<prec_t, 9> box;
    vector<prec_t> mass;
    vector<index_t> atype;
    vector<prec_t> coords;
    vector<bool_t> moves;
    try
    {
        if (!in_stru)
            throw "Can't find the structure file.";
        in_stru.clear();
        in_stru.seekg(0); // back to position 0
        
        if (in_stru.good())
        {
            char c_header[129];
            in_stru.getline(c_header, 129, '\n');
            header = c_header;
        }
        else
            throw "Invalid header.";

        if (in_stru.good())
            in_stru >> scaling;

        if (in_stru.good())
        {
            for (int i = 0 ; i < 9 ; ++i)
                in_stru >> box[i];
            in_stru.ignore(36, '\n');
        }

        if (in_stru.good())
        {
            vector<string> temp_elements;
            bool check = false;
            if (elements.size())
            {
                temp_elements.swap(elements);
                check = true;
            }
            string temp;
	    string line;
            getline(in_stru, line);
            istringstream iss(line);
            while (iss >> temp)
            {
                elements.push_back(temp);
                mass.push_back(ele_mass.at(temp));
            }
            /*
            while (in_stru.peek() != '\n')
            {
                if (in_stru.peek() == ' ' || in_stru.peek() == '\t')
                    in_stru.get();
                else
                {
                    in_stru >> temp;
                    elements.push_back(temp);
                    mass.push_back(ele_mass.at(temp));
                }
            }
            in_stru.get();
            if (check)
            {
                if (temp_elements.size() != elements.size())
                    throw "The numbers of elements are different is two structures.";
                for (int i = 0 ; i < elements.size() ; ++i)
                {
                    if (temp_elements[i] != elements[i])
                        throw "Different element order.";
                }
            }
	    */
        }

        if (in_stru.good())
        {
            size_t temp;
            index_t index = 0;
	    string line;
            getline(in_stru, line);
            istringstream iss(line);
            while (iss >> temp)
            {
                vector<index_t> temp_atype(temp, index);
                atype.insert(atype.end(), temp_atype.begin(), temp_atype.end());
                ++index;
            }
            /*
            while (in_stru.peek() != '\n')
            {
                if (in_stru.peek() == ' ' || in_stru.peek() == '\t')
                    in_stru.get();
                else
                {
                    in_stru >> temp;
                    vector<index_t> temp_atype(temp, index);
                    atype.insert(atype.end(), temp_atype.begin(), temp_atype.end());
                    ++index;
                }
            }
            in_stru.get();
	    */
        }

        move_type m_type;
        pos_type p_type;
        if (in_stru.good())
        {
            char coord_type[36];
            in_stru.getline(coord_type, 36, '\n');
            if (coord_type[0] == 'S')
            {
                m_type = SELECTIVE;
                in_stru.getline(coord_type, 36, '\n');
            }
            else
                m_type = FREE;
            if (coord_type[0] == 'C')
                p_type = CARTESIAN;
            else if (coord_type[0] == 'D')
                p_type = DIRECT;
            else
                throw "Invalid position coordinates type.";
        }

        if (in_stru.good())
        {
            prec_t x, y, z;
            char m_x, m_y, m_z;
            vector<prec_t> temp_coords;
            vector<bool_t> temp_moves(3, true);
            for (int i = 0 ; i < atype.size() ; ++i)
            {
                in_stru >> x >> y >> z;
                if (p_type == CARTESIAN)
                    temp_coords = {x, y, z};
                else if (p_type == DIRECT)
                    temp_coords = vector<prec_t>{x, y, z} * box;
                coords.insert(coords.end(), temp_coords.begin(), temp_coords.end());
                if (m_type == SELECTIVE)
                {
                    in_stru >> m_x >> m_y >> m_z;
                    temp_moves[0] = m_x == 'T' ? true : false;
                    temp_moves[1] = m_y == 'T' ? true : false;
                    temp_moves[2] = m_z == 'T' ? true : false;
                }
                moves.insert(moves.end(), temp_moves.begin(), temp_moves.end());
            }
        }

        cell in_cell(box, atype, mass, coords, moves, p_type, m_type);
        return in_cell;
    }
    catch (const char *msg)
    {
        cerr << msg << "\n";
    }
}

void KMC::in_out::print_stru(const cell &Cell, const string &stru_file, int precision) const
{
    ofstream out_stru(stru_file);
    out_stru.setf(ios::fixed);
    out_stru.setf(ios::showpoint);
    out_stru << header << "\n";
    out_stru << setprecision(2) << scaling << "\n";
    out_stru << setprecision(6);
    array<prec_t, 9> temp_box = Cell.get_box();
    double lx = temp_box[0],ly = temp_box[4],lz = temp_box[8];
    for (int i = 0 ; i < 3 ; ++i)
        out_stru << temp_box[i*3] << " " << temp_box[i*3+1] << " " << temp_box[i*3+2] << "\n";
    for (auto ele : elements)
        out_stru << ele << " ";
    out_stru << "\n";
    unordered_map<index_t, size_t> ele_num;
    for (int i = 0 ; i < elements.size() ; ++i)
        ele_num.insert(pair<index_t, size_t>(i, 0));
    for (auto type : Cell.get_atype())
        ele_num[type] += 1;
    for (int i = 0 ; i < elements.size() ; ++i)
        out_stru << ele_num[i] << " ";
    out_stru << "\n";
    out_stru << setprecision(precision);
    const atom_vec<bool_t> &moves = Cell.get_moves();
    const array<prec_t, 9> inv_box = inverse(temp_box);
    if (Cell.m_type == SELECTIVE)
        out_stru << "Selective\n";
    if (Cell.p_type == CARTESIAN)
        out_stru << "Cartesian\n";
    else if (Cell.p_type == DIRECT)
        out_stru << "Direct\n";
    vector<prec_t> temp_pos;
    for (int i = 0 ; i < Cell.get_atom_num() ; ++i)
    {
        if (Cell.p_type == CARTESIAN)
	{
	    temp_pos = vector<prec_t>{Cell.positions[i][0], Cell.positions[i][1], Cell.positions[i][2]};
	    if (temp_pos[0] < 0.0) temp_pos[0] = temp_pos[0] + lx; 
            if (temp_pos[1] < 0.0) temp_pos[1] = temp_pos[1] + ly; 
            if (temp_pos[2] < 0.0) temp_pos[2] = temp_pos[2] + lz; 
            if (temp_pos[0] > lx)  temp_pos[0] = temp_pos[0] - lx; 
            if (temp_pos[1] > ly)  temp_pos[1] = temp_pos[1] - ly; 
            if (temp_pos[2] > lz)  temp_pos[2] = temp_pos[2] - lz; 
            //out_stru << Cell.positions[i][0] << " " << Cell.positions[i][1] << " " << Cell.positions[i][2];
	    out_stru << temp_pos[0] << " " << temp_pos[1] << " " << temp_pos[2];
	}
        else if (Cell.p_type == DIRECT)
        { 
            temp_pos = vector<prec_t>{Cell.positions[i][0], Cell.positions[i][1], Cell.positions[i][2]} * inv_box;
            out_stru << temp_pos[0] << " " << temp_pos[1] << " " << temp_pos[2];
        }
        if (Cell.m_type == SELECTIVE)
            out_stru << " " << (moves[i][0] ? 'T' : 'F') << " " << (moves[i][1] ? 'T' : 'F') << " " << (moves[i][2] ? 'T' : 'F');
        out_stru << "\n";
    }
    out_stru.close();
}
