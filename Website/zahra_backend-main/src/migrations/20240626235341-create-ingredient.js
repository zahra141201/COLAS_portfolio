'use strict';
/** @type {import('sequelize-cli').Migration} */
module.exports = {
  async up(queryInterface, Sequelize) {
    await queryInterface.createTable('Ingredients', {
      id: {
        allowNull: false,
        autoIncrement: true,
        primaryKey: true,
        type: Sequelize.INTEGER
      },
      name: {
        type: Sequelize.STRING
      },
      expiration_date: {
        type: Sequelize.DATE
      },
      weight: {
        type: Sequelize.INTEGER
      },
      bought_date: {
        type: Sequelize.DATE
      },
      price: {
        type: Sequelize.FLOAT
      },
      description: {
        type: Sequelize.STRING
      },
      owner: {
        type: Sequelize.STRING,
        references: { model: 'Users', key: 'email' },
        onUpdate: 'CASCADE',
        onDelete: 'CASCADE',  // Ajout de la suppression en cascade
        allowNull: false,
      },
      createdAt: {
        allowNull: false,
        type: Sequelize.DATE
      },
      updatedAt: {
        allowNull: false,
        type: Sequelize.DATE
      }
    });
  },
  async down(queryInterface, Sequelize) {
    await queryInterface.dropTable('Ingredients');
  }
};
