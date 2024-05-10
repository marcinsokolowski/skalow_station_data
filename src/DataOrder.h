/*
 * DataOrder.h
 *
 *  Created on: 16 Apr 2018
 *      Author: ban115
 */

#ifndef DATAORDER_H_
#define DATAORDER_H_

// enum class DataOrder { TFBP, FTBP, BPFT, BPTF };
enum  DataOrder { DataOrder_TFBP=0, DataOrder_FTBP=1, DataOrder_BPFT=2, DataOrder_BPTF=3 };

DataOrder data_order_from_string(const char* str);

#endif /* DATAORDER_H_ */
