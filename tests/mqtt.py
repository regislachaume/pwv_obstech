from pwv_obstech.comm.mqtt import obstech_mqtt_client

client = obstech_mqtt_client(topic='Test')
client.publish(dict(test='test'))
