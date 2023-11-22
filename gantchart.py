import matplotlib.pyplot as plt


# Define tasks and their start and end dates
tasks = [
    {"task": "Task 1", "start": "2023-10-25", "end": "2023-11-05"},
    {"task": "Task 2", "start": "2023-11-06", "end": "2023-11-15"},
    {"task": "Task 3", "start": "2023-11-16", "end": "2023-11-30"},
]

# Create a figure and axis
fig, ax = plt.subplots(figsize=(10, 4))

# Plot the tasks as horizontal bars
for i, task in enumerate(tasks):
    start_date = task["start"]
    end_date = task["end"]
    ax.barh(task["task"], left=start_date, width=2, height=0.4, color=f"C{i}")

# Set the x-axis limits based on the task dates
start_date = min(task["start"] for task in tasks)
end_date = max(task["end"] for task in tasks)
ax.set_xlim(start_date, end_date)

# Customize the appearance
ax.set_xlabel("Timeline")
ax.set_title("Gantt Chart")
ax.invert_yaxis()  # Display tasks from top to bottom

# Show the Gantt chart
plt.tight_layout()
plt.show()