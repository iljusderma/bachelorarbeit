ilja@tss24:~$ ip address
1: lo: <LOOPBACK,UP,LOWER_UP> mtu 65536 qdisc noqueue state UNKNOWN group default qlen 1000
    link/loopback 00:00:00:00:00:00 brd 00:00:00:00:00:00
    inet 127.0.0.1/8 scope host lo
       valid_lft forever preferred_lft forever
    inet6 ::1/128 scope host 
       valid_lft forever preferred_lft forever
2: enp0s25: <BROADCAST,MULTICAST,UP,LOWER_UP> mtu 1500 qdisc fq_codel state UP group default qlen 1000
    link/ether 90:1b:0e:96:a1:c7 brd ff:ff:ff:ff:ff:ff
    inet 141.30.17.84/24 brd 141.30.17.255 scope global enp0s25
       valid_lft forever preferred_lft forever
    inet6 fe80::921b:eff:fe96:a1c7/64 scope link 
       valid_lft forever preferred_lft forever

ilja@tss24:~$ hostname -I
141.30.17.84 
