import serial
import csv
from datetime import datetime

# Configure the serial port
ser = serial.Serial('COM5', 115200)  # Replace 'COM3' with your Arduino's port
csv_file = open('Tests/HMotorTest4Actual.csv', 'w', newline='')
csv_writer = csv.writer(csv_file)

# Write the CSV header
csv_writer.writerow(['Timestamp (UTC)', 'trval', 'tqP', 'tq1P', 'tq2P'])

try:
    while True:
        line = ser.readline().decode('utf-8').strip()  # Read a line from the serial port
        data = line.split(',')  # Split the line into CSV fields

        now = datetime.now()
        timestamp = now.strftime('%H:%M:%S') + f".{now.microsecond //1000:03d}"
        data.insert(0, timestamp)
        csv_writer.writerow(data)  # Write the data to the CSV file
        print(','.join(data))  # Optional: Print to console for debugging
except KeyboardInterrupt:
    print("Logging stopped.")
finally:
    ser.close()
    csv_file.close()