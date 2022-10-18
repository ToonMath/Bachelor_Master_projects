import numpy as np
import statistics
import logging
import time
from optibook.synchronous_client import Exchange
from scipy.stats import norm
from datetime import datetime,date
from math import log,sqrt,pi,exp
logger = logging.getLogger('client')
logger.setLevel('ERROR')

print("Setup was successful.")

#connect to exchange
e = Exchange()
a = e.connect()
expire_time = "2021_04_16 12:00"
r = 0.0
sigma = 3.0

# Delete all outstanding orders
def Delete_all(instrument_id):
    outstanding = e.get_outstanding_orders(instrument_id)
    for o in outstanding.values():
        result = e.delete_order(instrument_id, order_id=o.order_id)
        print(f"Deleted order id {o.order_id}: {result}")


#calculate fair option price
def Fair_price(instrument_id):
    k=0
    total_book = []
    while len(e.get_last_price_book(instrument_id).bids)>0 and k<=1:
        book_bids = [0 for i in range(len(e.get_last_price_book(instrument_id).bids))]
        book_asks = [0 for i in range(len(e.get_last_price_book(instrument_id).asks))]
        j=0
        for j in range(len(e.get_last_price_book(instrument_id).bids)):
            book_bids[j] = e.get_last_price_book(instrument_id).bids[j].price
        for i in range(len(e.get_last_price_book(instrument_id).asks)):
            book_asks[i] = e.get_last_price_book(instrument_id).asks[i].price
        total_book = np.hstack((book_bids,book_asks))
        average = np.mean(total_book)
        zero_axis = average-statistics.median(total_book)
        dev = statistics.stdev(total_book)
      
        while dev>1 and np.abs(zero_axis)>1:
            if zero_axis>0:
                total_book = np.delete(total_book,np.argmin(total_book))
            elif zero_axis<0:
                total_book = np.delete(total_book,np.argmax(total_book))
            average = np.mean(total_book)
            zero_axis = -average+statistics.median(total_book)
            dev = statistics.stdev(total_book)
            k+=1
        return average
        
        #implement black scholes
def d1(S,K,T,r,sigma):
    return (log(S/K)+r*T+0.5*T*sigma**2)/(sigma*sqrt(T))
def d2(S,K,T,r,sigma):
    return d1(S,K,T,r,sigma)-sigma*sqrt(T)
def call_value(S,K,T,r,sigma):
    return S*norm.cdf(d1(S,K,T,r,sigma))-K*exp(-r*T)*norm.cdf(d2(S,K,T,r,sigma))

def put_value(S,K,T,r,sigma):
    return K*exp(-r*T)-S+call_value(S,K,T,r,sigma)

def calculate_time_to_date():
    today = datetime.now()
    t = (datetime.strptime(expire_time,'%Y_%m_%d %H:%M') - datetime.utcnow()).days / 365
    return t
def desired_bid(K,instrument_id):
    T = calculate_time_to_date()
    S = Fair_price(instrument_id)
    r = 0.0
    sigma = 3.0
    bid = call_value(S,K,T,r,sigma)+0.2
    return bid
def desired_ask(K,instrument_id):
    T = calculate_time_to_date()
    S = Fair_price(instrument_id)
    r = 0.0
    sigma = 3.0
    ask = put_value(S,K,T,r,sigma)+1.0
    return ask
    
def stock_list():
    instrumentlist = []
    positions = e.get_positions()
    for p in positions:
        instrumentlist.append(p)
    return instrumentlist
    
Strike_price = [50, 75, 100, 125,150,50, 75, 100, 125,150]

#determine if constraints are reached
def total_position(p):
    Total = 0
    positions = e.get_positions()
    Not_Max_positions = True
    if positions[p]>90 or positions[p]<-90:
        Max_positions = True
        
    for p in positions:
       Total+=positions[p]
    time.sleep(0.2)
    return Max_positions, positions[p],Total
    
def total_orders(p):
    Orders = 0
    Max_orders = False
    orders = e.get_outstanding_orders(p)
    for o in orders.values():
        Orders+=np.abs(o.volume)
    if Max_orders >=190:
        Max_orders = True
    return Max_orders
    
def call_delta(S,K):
    T = calculate_time_to_date()
    return norm.cdf(d1(S,K,T,r,sigma))
def put_delta(S,K):
    T = calculate_time_to_date()
    return -norm.cdf(-d1(S,K,T,r,sigma))

def Total_delta():
    Delta = 0
    substring_1 = "Call"
    substring_2 = "Put"
    for p in positions: 
            if substring_1 in p:
                Delta+=positions[p]*call_delta(Fair_price(p),Strike_price[labels.index(p)-1]) 
                time.sleep(0.1)
            elif substring_2 in p:
                Delta+=positions[p]*put_delta(Fair_price(p),Strike_price[labels.index(p)-1])
                time.sleep(0.1)
            else:
                Delta+= positions[p]
    print(Delta)
    return Delta
#start trading

cycle =1000
labels = stock_list()
for i in range(cycle):
    Vol = 5
    positions = e.get_positions()
    for p in positions:
        if p !="BMW" and not total_orders(p) and total_position(p)[1]+Vol<90:
            strike =Strike_price[labels.index(p)-1]
            buy = float(round(desired_bid(strike,p),2))
            result = e.insert_order(p, price=buy, volume=Vol, side='ask', order_type='limit')
            print(f"Order Id: {result}")
            trades = e.poll_new_trades(p)
            while not trades:
                print("No trades yet")
                time.sleep(0.2)
        elif p !="BMW" and not total_orders(p) and total_position(p)[1]-Vol>-90 and total_position(p)[2]>-500:
            strike =Strike_price[labels.index(p)-1]
            sell = float(round(desired_ask(strike,p),2))
            result = e.insert_order(p, price=buy, volume=Vol, side='bid', order_type='limit')
            print(f"Order Id: {result}")
            trades = e.poll_new_trades(p)
            while not trades:
                print("No trades yet")
                time.sleep(0.2)
        
            
        print("Now hedging")
    #perform the hedge
        Del = float(round(Total_delta()))
        print(Del)
        if Del<=-4 and not total_position("BMW")[1]<90 and not total_orders("BMW"):
            Buy_volume = -1*Del
            buy = float(round(Fair_price(p),2))
            if (total_position("BMW")[1]+Buy_volume)<90:
                result = e.insert_order("BMW", price=buy, volume=Buy_volume, side='ask', order_type='limit')
                print(result)
        elif Del>=4 and total_position("BMW")[1]>-90 and not total_orders("BMW"):
            Sell_volume = Del
            sell = float(round(Fair_price(p),2))
            if (total_position("BMW")[1]-Sell_volume)>-90:
                result = e.insert_order("BMW", price=sell, volume=Sell_volume, side='bid', order_type='limit')
                print(result)
        time.sleep(1)
# deleting all outstanding orders
    print("Deleting remaining orders")
    for p in labels:
        Delete_all(p)
        
#display executed trades
    for trade in labels:
        trades = e.poll_new_trades(trade)
        if trades:
            time.sleep(0.5)
            for t in trades:
                print(f"[TRADED {t.instrument_id}] price({t.price}), volume({t.volume}), side({t.side})")
#print current positions
    print('Cycle complete: Status')
    positions = e.get_positions()
    for p in positions:
        print(p, positions[p])
    time.sleep(1)  
    print(f'P/L:{e.get_pnl()}')
    print('End of cycle.')
    