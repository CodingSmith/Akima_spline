/*
 * akima_spline.h
 *
 * akima spline interpolation library without external
 * dependencies
 *
 * ---------------------------------------------------------------------
 * Copyright (C) 2018, CodingSmith (wangbfslingman at gmail.com)
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 *
 */


#ifndef AKIMA_SPLINE_H
#define AKIMA_SPLINE_H

#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <iterator>

// unnamed namespace only because the implementation is in this
// header file and we don't want to export symbols to the obj files
namespace
{

namespace cs
{

class AkimaSpline {
public:
    AkimaSpline(){};
    AkimaSpline(const std::vector<float> &x, const std::vector<float> &y){
                    set_points(x, y);
                };

    void set_points(const std::vector<float> &x, const std::vector<float> &y);
    bool generate_points(std::vector<float> &x, std::vector<float> &y);
    bool constructed_flag;
    
private:
    std::vector<float> raw_slopes, argmented_slopes;
    std::vector<float> m_x, m_y;
    std::vector<float> m_b, m_c, m_d;
    std::vector<int> num_count, bin_idx;
    std::vector<float> diff_x, diff_y;
    std::vector<float> diff_argmented_slopes;
    std::vector<float> f1, f2, f12;
    int data_length;

    void diff(std::vector<float> &val, std::vector<float> &diff_val);
    bool validity_check();
    void hist_count(std::vector<float> &x, std::vector<int> &num_count, std::vector<int> &bin, std::vector<int> &bb);
};



void AkimaSpline::set_points(const std::vector<float> &x, const std::vector<float> &y){  
    m_x = x;
    m_y = y;
    data_length = m_x.size();

    if(!validity_check()){
        constructed_flag = false;  
    } else {
        for(size_t i = 0; i < diff_x.size(); ++i){
            float slope_temp = diff_y[i]/diff_x[i];
            raw_slopes.push_back(slope_temp);
        }

        float mm = 2*raw_slopes[0] - raw_slopes[1];
        float mmm = 2*mm - raw_slopes[0];
        float mp = 2*raw_slopes.back() - raw_slopes[raw_slopes.size() - 2];
        float mpp = 2*mp - raw_slopes.back();

        argmented_slopes.push_back(mmm);
        argmented_slopes.push_back(mm);
        argmented_slopes.insert(argmented_slopes.end(), raw_slopes.begin(), raw_slopes.end());
        argmented_slopes.push_back(mp);
        argmented_slopes.push_back(mpp);

        
        diff(argmented_slopes, diff_argmented_slopes);
        for(size_t i = 0; i < diff_argmented_slopes.size(); ++i){
            diff_argmented_slopes[i] = fabs(diff_argmented_slopes[i]);
        }

        f1.assign(diff_argmented_slopes.begin()+2, diff_argmented_slopes.end());
        f2.assign(diff_argmented_slopes.begin(), diff_argmented_slopes.end()-2);
        for(size_t i = 0; i < f1.size(); ++i){
            float slope_temp = f1[i] + f2[i];
            f12.push_back(slope_temp);
        }

        std::vector<float>::iterator biggest = std::max_element(std::begin(f12), std::end(f12));  
        std::vector<int> b_idx;
        for(size_t i = 0; i < f12.size(); ++i){
            if(f12[i] > 1e-8 * (*biggest)){
                b_idx.push_back(i);
            }
        }

        for(size_t i = 1; i <= m_y.size(); ++i){
            m_b.push_back(argmented_slopes[i]);
        }

        for(size_t i = 0; i < b_idx.size(); ++i){
            m_b[b_idx[i]] = (f1[b_idx[i]]*argmented_slopes[b_idx[i]+1] + \
                             f2[b_idx[i]]*argmented_slopes[b_idx[i]+2]) \
                             / f12[b_idx[i]];
        }

        for(size_t i = 0; i < raw_slopes.size(); ++i){
            float m_c_temp;
            m_c_temp = (3*raw_slopes[i] - 2*m_b[i] - m_b[i+1]) / diff_x[i];
            m_c.push_back(m_c_temp);
        }

        for(size_t i = 0; i < raw_slopes.size(); ++i){
            float m_d_temp;
            m_d_temp = (m_b[i] + m_b[i+1]- 2*raw_slopes[i])/pow(diff_x[i],2);
            m_d.push_back(m_d_temp);
        }

        constructed_flag = true;
    }
}

bool AkimaSpline::validity_check(){
    if(m_x.size() != m_y.size()){
        std::cerr << "The sizes of input x and y are not equal!" << std::endl;
        return false;
    } else if(m_x.size() <3 ){
        std::cerr << "The sizes of input x and y vectors are too short!" << std::endl;
        return false;
    }

    diff(m_x, diff_x);
    for(size_t i = 0; i< diff_x.size(); ++i){
        if(diff_x[i] <= 0){
            std::cerr << "The input x vector must be in strictly ascending order!" << std::endl;
            return false; 
        }
    }
    diff(m_y, diff_y);
    return true;
}

void AkimaSpline::diff(std::vector<float> &val, std::vector<float> &diff_val){
    for(size_t i = 0; i < val.size()-1; ++i){
       float temp = val[i+1] - val[i];
       diff_val.push_back(temp);  
    }
}    
    
bool AkimaSpline::generate_points(std::vector<float> &x, std::vector<float> &y){
    
    for(size_t i = 0; i < x.size(); ++i){
        if(x[i]<m_x.front() || x[i]>m_x.back()){
            std::cerr << "All interpolation points xi must lie between m_x.front() and m_x.back()" << std::endl;
            return false;
        }
    }

    std::vector<int> num_count, bin, bb;
    hist_count(x, num_count, bin, bb);

    std::vector<float> m_w;
    for(size_t i = 0; i < x.size(); ++i){
        float wj_temp = x[i] - m_x[bb[i]];
        m_w.push_back(wj_temp);
    }

    for(size_t i = 0; i < x.size(); ++i){
        float y_temp = ((m_w[i] * m_d[bb[i]] + m_c[bb[i]])*m_w[i] + m_b[bb[i]]) * m_w[i] + m_y[bb[i]];
        y[i] = y_temp;
    }

    return true;
}

void AkimaSpline::hist_count(std::vector<float> &x, std::vector<int> &num_count, std::vector<int> &bin, std::vector<int> &bb){
    
    for(size_t i = 0; i < m_x.size(); ++i){
        num_count.push_back(0);
    }
    for(size_t i = 0; i < x.size(); ++i){
        float dist = 1000;
        int belong_id = 1000;
        for(size_t j = 0; j < m_x.size(); ++j){
            float dist_temp = fabs(x[i] - m_x[j]);
            if(dist_temp<dist){
                dist = dist_temp;
                belong_id = j;
            }
        }
        bin.push_back(belong_id);
        num_count[belong_id] = num_count[belong_id] + 1;

    }

    for(size_t i = 0; i < bin.size(); ++i){
        bin[i] = bin[i]>(num_count.size()-1) ? (num_count.size()-1) : bin[i]; 
    }


    for(size_t i = 0; i < x.size(); ++i){
        int bb_temp = bin[i];
        bb.push_back(bb_temp);
    }
}





} // namespace cs


} // namespace

#endif /* TK_SPLINE_H */
