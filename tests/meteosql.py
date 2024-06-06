from mysql import connector

config = dict(
    host='10.0.11.3',
    database='ElSauceWeather',
    user='obstech', 
    password='obstech4860'
)
cnx = connector.connect(**config)
cursor = cnx.cursor()

query = (
    "SELECT fecha, pressure, humidity, temperature_external "
    "FROM weather " 
    "WHERE fecha BETWEEN '2024-06-05T00:00:00' AND '2024-06-05T00:15:00' "
    "ORDER by fecha"
)
cursor.execute(query)
for (t, P, h, T) in cursor:
    print(f"{t} {P:6.1f} {T:+5.1f} {h:4.1f}")
