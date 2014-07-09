/**********************************************************************
 * Skewer - a fast and accurate adapter trimming tool
 *          using the bit-masked k-difference matching algorithm
 * Copyright (c) 2013-2014 by Hongshan Jiang
 * hongshan.jiang@gmail.com
 *
 * If you use this program, please cite the paper:
 * Jiang, H., Lei, R., Ding, S.W. and Zhu, S. (2014) Skewer: a fast and
 * accurate adapter trimmer for next-generation sequencing paired-end reads.
 * BMC Bioinformatics, 15, 182.
 * http://www.biomedcentral.com/1471-2105/15/182
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
#include <stdio.h>
#include <float.h>
#include <algorithm>
#include <string.h>
#include "matrix.h"
#include "fastq.h"

CODE map[256] = {
	//  0        1        2        3        4        5        6        7        8        9        A        B        C        D        E        F
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // 0
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // 1
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // 2
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // 3
	// 0x41~0x5A, A~Z
	CD_NONE,    CD_A,    CD_B,    CD_C,    CD_D, CD_NONE, CD_NONE,    CD_G,    CD_H, CD_NONE, CD_NONE,    CD_K, CD_NONE,    CD_M,    CD_N, CD_NONE, // 4
	CD_NONE, CD_NONE,    CD_R,    CD_S,    CD_T,    CD_T,    CD_V,    CD_W, CD_NONE,    CD_Y, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // 5
	// 0x61~0x7A, a~z
	CD_NONE,    CD_A,    CD_B,    CD_C,    CD_D, CD_NONE, CD_NONE,    CD_G,    CD_H, CD_NONE, CD_NONE,    CD_K, CD_NONE,    CD_M,    CD_N, CD_NONE, // 6
	CD_NONE, CD_NONE,    CD_R,    CD_S,    CD_T,    CD_T,    CD_V,    CD_W, CD_NONE,    CD_Y, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // 7

	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // 8
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // 9
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // A
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // B
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // C
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // D
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // E
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE  // F
};

bool blurry[256] = {
	// 0     1     2     3     4     5     6     7     8     9     A     B     C     D     E     F
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // 0
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // 1
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // 2
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // 3

	true,false, true,false, true, true, true,false, true, true, true, true, true, true, true, true, // 4
	true, true, true, true,false,false, true, true, true, true, true, true, true, true, true, true, // 5

	true,false, true,false, true, true, true,false, true, true, true, true, true, true, true, true, // 6
	true, true, true, true,false,false, true, true, true, true, true, true, true, true, true, true, // 7

	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // 8
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // 9
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // A
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // B
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // C
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // D
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // E
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true  // F
};

CODE complement[CD_CNT] = {
	CD_NONE, CD_T, CD_G, CD_C, CD_A, CD_Y, CD_R, CD_W, CD_S, CD_M, CD_K, CD_V, CD_H, CD_D, CD_B, CD_N
};

char character[CD_CNT] = {
	'N', 'A', 'C', 'G', 'T', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N'
};

double scoring[CD_CNT][CD_CNT] = {
	//   -    A    C    G    T  | R    Y    S    W    K    M |  B    D    H    V |  N
	{    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, 0.05}, //-
	{    1,   0,   1,   1,   1,   0,   1,   1,   0,   1,   0,   1,   0,   0,   0, 0.05}, //A
	{    1,   1,   0,   1,   1,   1,   0,   0,   1,   1,   0,   0,   1,   0,   0, 0.05}, //C
	{    1,   1,   1,   0,   1,   0,   1,   0,   1,   0,   1,   0,   0,   1,   0, 0.05}, //G
	{    1,   1,   1,   1,   0,   1,   0,   1,   0,   0,   1,   0,   0,   0,   1, 0.05}, //T

	{    1,   0,   1,   0,   1,   0,   1,0.75,0.75,0.75,0.75, 0.5,   0, 0.5,   0, 0.05}, //R
	{    1,   1,   0,   1,   0,   1,   0,0.75,0.75,0.75,0.75,   0, 0.5,   0, 0.5, 0.05}, //Y
	{    1,   1,   0,   0,   1,0.75,0.75,   0,   1,0.75,0.75,   0, 0.5, 0.5,   0, 0.05}, //S
	{    1,   0,   1,   1,   0,0.75,0.75,   1,   0,0.75,0.75, 0.5,   0,   0, 0.5, 0.05}, //W
	{    1,   1,   1,   0,   0,0.75,0.75,0.75,0.75,   0,   1,   0,   0, 0.5, 0.5, 0.05}, //K
	{    1,   0,   0,   1,   1,0.75,0.75,0.75,0.75,   1,   0, 0.5, 0.5,   0,   0, 0.05}, //M

	{    1,   1,   0,   0,   0, 0.5,   0,   0, 0.5,   0, 0.5,   0, 0.4, 0.4, 0.4, 0.05}, //B
	{    1,   0,   1,   0,   0,   0, 0.5, 0.5,   0,   0, 0.5, 0.4,   0, 0.4, 0.4, 0.05}, //D
	{    1,   0,   0,   1,   0, 0.5,   0, 0.5,   0, 0.5,   0, 0.4, 0.4,   0, 0.4, 0.05}, //H
	{    1,   0,   0,   0,   1,   0, 0.5,   0, 0.5, 0.5,   0, 0.4, 0.4, 0.4,   0, 0.05}, //V

	{ 0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,    0}  //N
};

uint64 chrVadp[CD_CNT][CD_CNT] = {
// adp   -    A    C    G    T  | R    Y    S    W    K    M |  B    D    H    V |  N   //chr
	{    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0}, //-
	{    1,   0,   1,   1,   1,   0,   1,   1,   0,   1,   0,   1,   0,   0,   0,   0}, //A
	{    1,   1,   0,   1,   1,   1,   0,   0,   1,   1,   0,   0,   1,   0,   0,   0}, //C
	{    1,   1,   1,   0,   1,   0,   1,   0,   1,   0,   1,   0,   0,   1,   0,   0}, //G
	{    1,   1,   1,   1,   0,   1,   0,   1,   0,   0,   1,   0,   0,   0,   1,   0}, //T

	{    1,   1,   1,   1,   1,   0,   1,   1,   1,   1,   1,   1,   0,   1,   0,   0}, //R
	{    1,   1,   1,   1,   1,   1,   0,   1,   1,   1,   1,   0,   1,   0,   1,   0}, //Y
	{    1,   1,   1,   1,   1,   1,   1,   0,   1,   1,   1,   0,   1,   1,   0,   0}, //S
	{    1,   1,   1,   1,   1,   1,   1,   1,   0,   1,   1,   1,   0,   0,   1,   0}, //W
	{    1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   1,   0,   0,   1,   1,   0}, //K
	{    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   1,   1,   0,   0,   0}, //M

	{    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   1,   1,   1,   0}, //B
	{    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   1,   1,   0}, //D
	{    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   1,   0}, //H
	{    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   0}, //V

	{    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0}  //N
};

const double MIN_PENALTY = 0.477121255;
const double MEAN_PENALTY = 2.477121255;
const double MAX_PENALTY = 4.477121255;
///////////////////////////////////////

bool cElementSet::insert (const ELEMENT& val)
{
	pair<ELEMENT_SET::iterator,bool> ret;
	ELEMENT_SET::iterator it = this->find(val);
	if(it == this->end()){
		ret = ELEMENT_SET::insert(val);
		return ret.second;
	}
	if(val.score < it->score){
		return true;
	}
	this->erase(it++);
	ELEMENT_SET::insert(it, val);
	return true;
}

///////////////////////////////////////
cAdapter::cAdapter()
{
	len = 0;
}

cAdapter::~cAdapter()
{
}

void cAdapter::Init(char * seq, size_t sLen, TRIM_MODE trimMode)
{
	int i;
	// construct sequence
	this->len = (int(sLen) > MAX_ADAPTER_LEN) ? MAX_ADAPTER_LEN : sLen;
	gzstrncpy(sequence, seq, len);
	this->trimMode = trimMode;

	// construct mismatch bits
	int code, code2;
	uint64 bits;
	for(code=0; code<CD_CNT; code++){
		bits = 0;
		for(i=int(len)-1; i>=0; i--){
			code2 = map[uchar(sequence[i])];
			bits = (bits << 1) | chrVadp[code][code2];
		}
		matchBits[code] = ~bits;
	}
}

void cAdapter::Init2(char * seq, size_t sLen)
{
	int i;
	// construct sequence
	if(int(sLen) > MAX_ADAPTER_LEN){
		seq += (sLen - MAX_ADAPTER_LEN);
		this->len = MAX_ADAPTER_LEN;
	}
	else{
		this->len = sLen;
	}
	for(i=0; i<int(len); i++)
		sequence[i] = character[complement[map[uchar(seq[len-1-i])]]];
	sequence[len] = '\0';
	this->trimMode = TRIM_TAIL;

	// construct mismatch bits
	int code, code2;
	uint64 bits;
	for(code=CD_BASIC_CNT-1; code>=0; code--){
		bits = 0;
		for(i=int(len)-1; i>=0; i--){
			code2 = map[uchar(sequence[i])];
			bits = (bits << 1) | chrVadp[code][code2];
		}
		matchBits[code] = ~bits;
	}
	for(bits=~bits,code=CD_BASIC_CNT; code<CD_CNT; code++){
		matchBits[code] = bits;
	}
}

inline void cAdapter::UPDATE_COLUMN(deque<ELEMENT> & queue, uint64 &d0bits, uint64 &lbits, uint64 &unbits, uint64 &dnbits, double &penal, double &dMaxPenalty, int &iMaxIndel)
{
	int i;
	double score;
	uint64 bits = ~lbits | d0bits;
	for(bits>>=1,i=1; i<int(queue.size())-1; i++,bits>>=1){
		if((bits & 0x01) == 0){
			if(cMatrix::bSensitive){
				score = queue[i].score + (penal - cMatrix::dDelta);
				if( (queue[i-1].score < score) && (queue[i-1].nIndel < iMaxIndel) ){
					if( (queue[i+1].score < score) && (queue[i+1].nIndel < iMaxIndel) ){
						if(queue[i-1].score < queue[i+1].score){
							queue[i] = queue[i-1];
							dnbits |= (1L << (i-1));
						}
						else{
							queue[i] = queue[i+1];
							unbits |= (1L << (i+1));
						}
					}
					else{
						queue[i] = queue[i-1];
						dnbits |= (1L << (i-1));
					}
					queue[i].nIndel++;
				}
				else{
					if( (queue[i+1].score < score) && (queue[i+1].nIndel < iMaxIndel) ){
						queue[i] = queue[i+1];
						unbits |= (1L << (i+1));
						queue[i].nIndel++;
					}
					else{
						queue[i].score = score;
					}
				}
				queue[i].score += cMatrix::dDelta;
			}
			else{ // !cMatrix::bSensitive
				queue[i].score += penal;
			}
			if(queue[i].score >= dMaxPenalty){
				lbits &= ~(1L << i);
			}
		}
	}
	if(queue.size() > 1){
		if((bits & 0x01) == 0){
			if(cMatrix::bSensitive){
				if( (queue[i-1].nIndel < iMaxIndel) && (queue[i-1].score + cMatrix::dDelta < queue[i].score + penal) ){
					queue[i] = queue[i-1];
					dnbits |= (1L << (i-1));
					queue[i].score += cMatrix::dDelta;
					queue[i].nIndel++;
				}
				else{
					queue[i].score += penal;
				}
			}
			else{ // !cMatrix::bSensitive
				queue[i].score += penal;
			}
			if(queue[i].score >= dMaxPenalty){
				lbits &= ~(1L << i);
			}
		}
		for(; i>0; i--){
			if(queue.back().score < dMaxPenalty) break;
			queue.pop_back();
		}
	}
}

bool cAdapter::align(char * read, size_t rLen, uchar * qual, size_t qLen, cElementSet &result, int bc, bool bBestAlign)
{
	bool bDetermined = false;
	ELEMENT elem;
	double dMaxPenalty = cMatrix::dPenaltyPerErr * len + 0.001;
	int minK = bBestAlign ? cMatrix::iMinOverlap : 1;
	int iMaxIndel = ceil(cMatrix::dEpsilonIndel * len);
	double dMu = (bc >= 0) ? cMatrix::dMu : MIN_PENALTY;

	deque<ELEMENT> queue;
	ELEMENT element;
	double score;
	uint64 legalBits = 0;
	int i, j, jj;
	element.idx.bc = bc;
	if(trimMode & TRIM_HEAD){
		for(i=1; i<=int(len)-minK; i++){
			element.idx.pos = -i;
			element.score = cMatrix::dPenaltyPerErr * i;
			element.nIndel = 0;
			queue.push_back(element);
			legalBits = (legalBits << 1) | 1;
		}
	}
	else{
		for(i=1,score=cMatrix::dDelta; i<int(len); i++,score+=cMatrix::dDelta){
			if(i > iMaxIndel) break;
			element.idx.pos = -i;
			element.score = score;
			element.nIndel = i;
			queue.push_back(element);
			legalBits = (legalBits << 1) | 1;
		}
	}
	element.nIndel = 0;
	uint64 mbits, xbits, unbits, dnbits, d0bits;
	unbits = dnbits = 0L;
	double penal;
	for(j=0; j<int(rLen); j++){
		jj = j;
		mbits = matchBits[map[uchar(read[jj])]];
		penal = ((qLen > 0) ? cMatrix::penalty[qual[jj]] : dMu);

		element.idx.pos = j;
		element.score = ((mbits & 0x01) == 0) ? penal : 0;
		queue.push_front(element);

		xbits = mbits | unbits;
		dnbits <<= 1;
		unbits <<= 1;
		d0bits = ((dnbits + (xbits & dnbits)) ^ dnbits) | xbits;
		legalBits = (legalBits << 1) | 1;

		UPDATE_COLUMN(queue, d0bits, legalBits, unbits, dnbits, penal, dMaxPenalty, iMaxIndel);

		dnbits &= d0bits;
		unbits &= d0bits;

		if(queue.size() == len){
			if(bBestAlign){
				if(trimMode == TRIM_HEAD){
					i = (queue.back().idx.pos < 0) ? (len + queue.back().idx.pos) : len;
					if( !bDetermined || (i * cMatrix::dMu - queue.back().score) > elem.score * (i+1) ){
						elem = queue.back();
						dMaxPenalty = elem.score;
						elem.score = (i * cMatrix::dMu - elem.score) / (i+1); // normalization
						elem.idx.pos = rLen - 1 - j;
					}
				}
				else{
					elem = queue.back();
					dMaxPenalty = elem.score;
					elem.score = (len * cMatrix::dMu - elem.score) / (len + 1); // normalization
				}
				bDetermined = true;
				if(dMaxPenalty == 0) break;
			}
			else{
				elem = queue.back();
				elem.score = len * cMatrix::dMu - elem.score; // normalization
				result.insert(elem);
			}
			queue.pop_back();
		}
	}
	if(dMaxPenalty > 0){ // not the case of "perfect match for single-end reads trimming"
		if(bBestAlign){
			if(trimMode & TRIM_TAIL){
				dMaxPenalty = (cMatrix::dPenaltyPerErr * queue.size() + 0.001);
				for(i=queue.size(); i>=minK; i--, dMaxPenalty-=cMatrix::dPenaltyPerErr){
					if(dMaxPenalty <= 0) break;
					if(queue.back().score < dMaxPenalty){
						if(!bDetermined || ((i * cMatrix::dMu - queue.back().score) > elem.score * (i+1)) ){
							elem = queue.back();
							dMaxPenalty = elem.score;
							elem.score = (i * cMatrix::dMu - elem.score) / (i+1); // normalization
							bDetermined = true;
						}
					}
					queue.pop_back();
				}
			}
			else{
				dMaxPenalty -= (len - queue.size()) * cMatrix::dDelta;
				iMaxIndel -= (len - queue.size());
				for(i=queue.size(); i>=minK; i--, dMaxPenalty-=cMatrix::dDelta, iMaxIndel--){
					if( (dMaxPenalty <= 0) || (iMaxIndel < 0) ) break;
					if( (queue.back().score < dMaxPenalty) && (queue.back().nIndel <= iMaxIndel) ){
						if(!bDetermined || ((i * cMatrix::dMu - queue.back().score) > elem.score * (i+1)) ){
							elem = queue.back();
							dMaxPenalty = elem.score;
							elem.score = (i * cMatrix::dMu - elem.score) / (i+1); // normalization
							elem.idx.pos = -(len - i);
							bDetermined = true;
						}
					}
					queue.pop_back();
				}
			}
		}
		else{
			dMaxPenalty = cMatrix::dPenaltyPerErr * queue.size() + 0.001;
			for(i=queue.size(); i>=minK; i--, dMaxPenalty-=cMatrix::dPenaltyPerErr){
				if(queue.back().score < dMaxPenalty){
					elem = queue.back();
					elem.score = i * cMatrix::dMu - elem.score; // normalization
					result.insert(elem);
				}
				queue.pop_back();
			}
		}
	}
	if(bDetermined){
		result.clear();
		result.insert(elem);
	}

	return bDetermined;
}

deque<cAdapter> cMatrix::firstAdapters;
deque<cAdapter> cMatrix::secondAdapters;
deque<cAdapter> cMatrix::junctionAdapters;
vector<int> cMatrix::junctionLengths;
bool cMatrix::bShareAdapter = false;
double cMatrix::dEpsilon = 0.15;
double cMatrix::dEpsilonIndel = 0.03;
double cMatrix::dPenaltyPerErr = cMatrix::dEpsilon * MEAN_PENALTY;
double cMatrix::dDelta = MAX_PENALTY;
double cMatrix::dMu = MEAN_PENALTY;
double cMatrix::penalty[256];
bool cMatrix::bSensitive = false;
int cMatrix::iMinOverlap = 3;

///////////////////////////////////////
cMatrix::cMatrix()
{
}

cMatrix::~cMatrix()
{
}

bool cMatrix::CalcRevCompScore(char * seq, char * seq2, int len, uchar * qual, uchar * qual2, size_t qLen, double &score)
{
	double dMaxPenalty = dPenaltyPerErr * len;
	double penal;
	CODE code, code2;
	score = 0.0;
	for(int i=0; i<len; i++){
		code = map[uchar(seq[i])];
		code2 = complement[map[uchar(seq2[len-1-i])]];
		penal = scoring[code][code2];
		if(penal > 0.0){
			if(qLen > 0){
				if(cMatrix::penalty[qual[i]] <= cMatrix::penalty[qual2[len-1-i]]){
					penal *= cMatrix::penalty[qual[i]];
				}
				else{
					penal *= cMatrix::penalty[qual2[len-1-i]];
				}
			}
			else{
				penal *= dMu;
			}
			score += penal;
			if(score > dMaxPenalty){
				return false;
			}
		}
	}
	score = len * dMu - score; // normalization
	return true;
}

//// public functions
void cMatrix::InitParameters(double dEpsilon, double dEpsilonIndel, int baseQual, bool bShareAdapter)
{
	cMatrix::dEpsilon = dEpsilon;
	cMatrix::dEpsilonIndel = dEpsilonIndel;
	cMatrix::dPenaltyPerErr = dEpsilon * MEAN_PENALTY;
	cMatrix::bSensitive = (dEpsilonIndel > 0);
	// pre-calcualte the penalties corresponding to quality values
	int chr;
	for(chr=0; chr<=baseQual; chr++){
		cMatrix::penalty[chr] = MIN_PENALTY;
	}
	int i;
	for(i=1; i<40; i++,chr++){
		cMatrix::penalty[chr] = MIN_PENALTY + i / 10.0;
	}
	for(; chr<256; chr++){
		cMatrix::penalty[chr] = MAX_PENALTY;
	}
	cMatrix::bShareAdapter = bShareAdapter;
	cMatrix::firstAdapters.clear();
	cMatrix::secondAdapters.clear();
	cMatrix::junctionAdapters.clear();
}

void cMatrix::AddAdapter(deque<cAdapter> & adapters, char * vector, size_t len, TRIM_MODE trimMode)
{
	cAdapter adapter;
	adapter.Init(vector, len, trimMode);
	adapters.push_back(adapter);
}

void cMatrix::CalculateJunctionLengths()
{
	deque<cAdapter>::iterator it_adapter;
	for(it_adapter=junctionAdapters.begin(); it_adapter!=junctionAdapters.end(); it_adapter++){
		junctionLengths.push_back((*it_adapter).len);
	}
}

bool cMatrix::isBlurry(char * seq, size_t len)
{
	size_t u;
	int iMaxBlurry = ceil(cMatrix::dEpsilon * len);
	int iBlurry = 0;
	for(u=0; u<len; u++){
		if(blurry[int(seq[u])]){
			if(++iBlurry > iMaxBlurry){
				return true;
			}
		}
	}
	return false;
}

bool cMatrix::checkQualities(uchar * quals, size_t len, int minQual)
{
	size_t u;
	if(len == 0) return true;
	int total = 0;
	for(u=0; u<len; u++){
		total += quals[u];
	}   
	return (double(total) / len) >= minQual;
}

int cMatrix::trimByQuality(uchar * quals, size_t len, int minQual)
{
	int i;
	for(i=(int)len-1; i>=0; i--){
		if(quals[i] >= minQual)
			break;
	}
	return (i+1);
}

INDEX cMatrix::findAdapter(char * read, size_t rLen, uchar * qual, size_t qLen)
{
	deque<cAdapter>::iterator it_adapter;
	cAdapter * pAdapter;
	cElementSet result;
	double maxScore = -1;
	INDEX index;
	index.pos = int(rLen);
	index.bc = -1;
	int i;
	for(i=0,it_adapter=firstAdapters.begin(); it_adapter!=firstAdapters.end(); it_adapter++,i++){
		pAdapter = &(*it_adapter);
		if(pAdapter->align(read, rLen, qual, qLen, result, i)){
			if(result.begin()->score > maxScore){
				index = result.begin()->idx;
				maxScore = result.begin()->score;
			}
		}
	}
	return index;
}

INDEX cMatrix::findJuncAdapter(char * read, size_t rLen, uchar * qual, size_t qLen)
{
	deque<cAdapter>::iterator it_adapter;
	cAdapter * pAdapter;
	cElementSet result;
	double maxScore = -1;
	INDEX index;
	index.pos = int(rLen);
	index.bc = -1;
	int i;
	for(i=0,it_adapter=junctionAdapters.begin(); it_adapter!=junctionAdapters.end(); it_adapter++,i++){
		pAdapter = &(*it_adapter);
		if(pAdapter->align(read, rLen, qual, qLen, result, i)){
			if(result.begin()->score > maxScore){
				index = result.begin()->idx;
				maxScore = result.begin()->score;
			}
		}
	}
	return index;
}

INDEX cMatrix::findAdapterWithPE(char * read, char * read2, size_t rLen, uchar * qual, uchar * qual2, size_t qLen)
{
	deque<cAdapter>::iterator it_adapter;
	cAdapter * pAdapter;
	cElementSet result;
	cElementSet result2;
	INDEX index;
	index.pos = int(rLen);
	index.bc = -1;
	int i;
	for(i=0,it_adapter=firstAdapters.begin(); it_adapter!=firstAdapters.end(); it_adapter++,i++){
		pAdapter = &(*it_adapter);
		pAdapter->align(read, rLen, qual, qLen, result, i, false);
	}
	deque<cAdapter> *pAdapters = (bShareAdapter ? &firstAdapters : &secondAdapters);
	for(i=0,it_adapter=pAdapters->begin(); it_adapter!=pAdapters->end(); it_adapter++,i++){
		pAdapter = &(*it_adapter);
		pAdapter->align(read2, rLen, qual2, qLen, result2, i, false);
	}
	if(result.empty()){
		if(result2.empty())
			return index;
		if(result2.begin()->idx.pos <= 0)
			return result2.begin()->idx;
	}
	else{
		if(result.begin()->idx.pos <= 0){
			index = result.begin()->idx;
			index.bc *= pAdapters->size();
			if(!result2.empty() && result2.begin()->idx.pos <= 0){
				index.bc += result2.begin()->idx.bc;
			}
			return index;
		}
	}
	double maxScore = -1;
	double score;
	cElementSet::iterator it_element, it_element2;
	bool bRevComplement;
	int pos, pos2, cpos;
	it_element = result.begin();
	it_element2 = result2.begin();
	while( true ){
		pos = (it_element == result.end()) ? INT_MAX : it_element->idx.pos;
		pos2 = (it_element2 == result2.end()) ? INT_MAX : it_element2->idx.pos;
		if( (pos == INT_MAX) && (pos2 == INT_MAX) )
			break;
		cpos = (pos <= pos2) ? pos : pos2;
		bRevComplement = CalcRevCompScore(read, read2, cpos, qual, qual2, qLen, score);
		if(pos < pos2){
			if(bRevComplement){
				score += it_element->score;
				if(score > maxScore){
					maxScore = score;
					index = it_element->idx;
					index.bc *= pAdapters->size();
				}
			}
			it_element++;
		}
		else if(pos > pos2){
			if(bRevComplement){
				score += it_element2->score;
				if(score > maxScore){
					maxScore = score;
					index = it_element2->idx;
				}
			}
			it_element2++;
		}
		else{ // ==
			if(bRevComplement){
				score += it_element->score + it_element2->score;
				if(score > maxScore){
					maxScore = score;
					index = it_element->idx;
					index.bc *= pAdapters->size();
					index.bc += it_element2->idx.bc;
				}
			}
			it_element++;
			it_element2++;
		}
	}

	return index;
}

INDEX cMatrix::mergePE(char * read, char * read2, size_t rLen, uchar * qual, uchar * qual2, size_t qLen, size_t startPos, size_t jLen)
{
	INDEX index;
	cElementSet result;
	cAdapter adapter;
	double score;
	int pos, clen;
	bool bRevComplement = false;
	size_t endPos, eLen;
	char chr;
	uchar uchr;
	int i;
	index.pos = int(rLen);
	index.bc = -1;
	adapter.Init2(read2, rLen);
	if(adapter.align(read, rLen, NULL, 0, result, -1)){
		pos = result.begin()->idx.pos;
		if(pos >= 0){
			clen = rLen - pos;
			bRevComplement = CalcRevCompScore(read + pos, read2 + pos, clen, qual + pos, qual2 + pos, qLen, score);
		}
	}
	if(bRevComplement){ // overlap detected
		if(qLen > 0)
			combinePairSeqs(read+pos, read2+pos, clen, qual+pos, qual2+pos, qLen);
		endPos = startPos + jLen;
		if(endPos + clen >= rLen){ // junction adapter locates in overlapping region
			index.pos = rLen - (endPos + clen - rLen);
		}
		else{
			eLen = rLen - (endPos + clen);
			index.pos += eLen;
			for(i=eLen/2; i>=0; i--){ //reverse
				chr = read2[endPos + i];
				read2[endPos + i] = read2[endPos + eLen - 1 - i];
				read2[endPos + eLen - 1 - i] = chr;
			}
			for(i=0; i<(int)eLen; i++){
				read2[startPos + i] = character[complement[map[uchar(read2[endPos + i])]]];
			}
			if(qLen > 0){
				for(i=eLen/2; i>=0; i--){ //reverse
					uchr = qual2[endPos + i];
					qual2[endPos + i] = qual2[endPos + eLen - 1 -i];
					qual2[endPos + eLen - 1 - i] = uchr;
				}
				for(i=0; i<(int)eLen; i++){
					qual2[startPos + i] = qual2[endPos + i];
				}
			}
		}
	}
	else{
		index.pos -= cMatrix::iMinOverlap;
	}
	return index;
}

void cMatrix::combinePairSeqs(char * read, char * read2, int len, uchar * qual, uchar * qual2, size_t qLen)
{
	CODE code, code2;
	for(int i=0; i<len; i++){
		code = map[uchar(read[i])];
		code2 = complement[map[uchar(read2[len-1-i])]];
		if(qual2[len-1-i] > qual[i]){
			qual[i] = qual2[len-1-i];
			if(code != code2){
				read[i] = character[code2];
			}
		}
	}
}