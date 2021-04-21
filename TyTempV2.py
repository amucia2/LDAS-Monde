import pandas as pd
import matplotlib.pyplot as plt


log = open('fahtemp.log','r')
df = pd.read_csv(log)
df.columns = ["date","temp","target","status"]
df['date'] = pd.to_datetime(df['date'],unit='s')

###df['temp']=df['temp'].rolling(21,min_periods=1).mean()
#temp7=pd.DataFrame(df['temp'].rolling(21,min_periods=1).mean())

ts = [0]
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylim(60,90)
i=0
for n in range(len(df['status'])-1):
    if df['status'].loc[n+1]==" FAH is running" and df['status'].loc[n]==" FAH is not running":
        ts.append(n)
        x = [df['temp'][ts[i]],df['temp'][n]]
        y = [df['date'][ts[i]],df['date'][n]]
        i+=1
        plt.plot(y,x,c='r')
        
    elif df['status'].loc[n+1]==" FAH is not running" and df['status'].loc[n]==" FAH is running":
        ts.append(n)
        x = [df['temp'][ts[i]],df['temp'][n]]
        y = [df['date'][ts[i]],df['date'][n]]
        i+=1
        plt.plot(y,x,c='b')

plt.plot(df['date'],df['target'],color='black')
plt.plot(df['date'],df['temp'].rolling(21,min_periods=1).mean(),color='cyan')
plt.tight_layout()     
plt.show() 