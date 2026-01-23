#include <stdio.h>
#include <stdlib.h>
#include <bsd/sys/queue.h>

// Structure for list entries containing integers
struct entry {
    int value;
    TAILQ_ENTRY(entry) entries;  // Tail queue pointers
};

// Define the head of the tail queue
TAILQ_HEAD(tailq_head, entry);

int main() {
    struct tailq_head head;
    struct entry *item, *temp;

    // Initialize the tail queue
    TAILQ_INIT(&head);

    printf("Creating TAILQ list with numbers 0-10...\n");

    // Create entries with values 0 to 10
    for (int i = 0; i <= 10; i++) {
        item = (struct entry *)malloc(sizeof(struct entry));
        if (item == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            return 1;
        }
        item->value = i;
        TAILQ_INSERT_TAIL(&head, item, entries);
        printf("Inserted: %d\n", i);
    }

    printf("\nList before removing odd numbers:\n");
    TAILQ_FOREACH(item, &head, entries) {
        printf("%d ", item->value);
    }
    printf("\n");

    printf("\nRemoving odd numbers using TAILQ_FOREACH_SAFE...\n");

    // Use TAILQ_FOREACH_SAFE to safely delete odd numbers
    TAILQ_FOREACH_SAFE(item, &head, entries, temp) {
        if (item->value % 2 != 0) {  // Check if odd
            printf("Removing: %d\n", item->value);
            TAILQ_REMOVE(&head, item, entries);
            free(item);
        }
    }

    printf("\nList after removing odd numbers (only even numbers remain):\n");
    TAILQ_FOREACH(item, &head, entries) {
        printf("%d ", item->value);
    }
    printf("\n");

    // Clean up: free remaining entries
    printf("\nCleaning up remaining entries...\n");
    while (!TAILQ_EMPTY(&head)) {
        item = TAILQ_FIRST(&head);
        TAILQ_REMOVE(&head, item, entries);
        free(item);
    }

    printf("Done!\n");

    return 0;
}
