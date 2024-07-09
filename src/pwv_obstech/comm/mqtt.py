from paho.mqtt import client as mqtt
from paho.mqtt.client import MQTTMessageInfo
import logging
from dataclasses import dataclass, field
import json
from typing import Union

@dataclass(frozen=True)
class MQTTClientInfo:
    host: str
    port: int
    username: str
    password: str
    topic: str 
    subscriptions: list[str] = field(default_factory=lambda: [])
    timeout: int = 60

class MQTTClient(mqtt.Client):

    def __init__(self, info: MQTTClientInfo) -> None: 
       
        self._info = info

        super().__init__(mqtt.CallbackAPIVersion.VERSION2)
        self.username_pw_set(info.username, password=info.password)

    def on_connect(self, userdata, flags, rc: int):

        host = self._info.host
        port = self._info.port 

        name = f"MQTT broker {host} on port {port}"
        if rc == 0:
            logging.info(f"Connected to {name}")
        else:
            msg = f"Failed to connect to {name}"
            logging.info(msg)
            raise ConnectionError(msg)

        for s in self._info.subscriptions:
            self.subscribe(s)

    def connect(self) -> int:

        host = self._info.host
        port = self._info.port
        timeout = self._info.timeout

        return super().connect(host, port, timeout)

    def publish(
        self, 
        payload: dict[str, object], 
        **kwargs
    ) -> MQTTMessageInfo:

        payload = json.dumps(payload, indent=4)
        topic = f"/ElSauce/Weather/{self._info.topic}"

        return super().publish(topic, payload, **kwargs)

def obstech_mqtt_client(
    topic,
    subscriptions: list[str] = []
) -> MQTTClient:

    info = MQTTClientInfo(
        host="10.0.11.3", 
        port=1883, 
        timeout=60, 
        username="obstech",
        password="obstech4860",
        subscriptions=subscriptions,
        topic=topic
    )
    return MQTTClient(info)
